"""
This module represents an explicit, direct sparse 2D finite-difference Poisson solver for heterogeneous media,
i.e. media having spatially varying (space-dependent) medium parameters.
The following functionality is currently available in this module:
* 'Forward' problem: given the source, and the conductivity distribution, solve for the electrical potential ψ:
	 Example problem: solve Aψ=∇⋅j; A=∇⋅(σ(x,z)∇), for ψ.
* 'Inverse' problem: given wave propagation related pore pressure, and mechanical medium properties, calculate the source term:  
	 Example problem: solve Aψ=-∇⋅j; A=∇⋅([Q*k/η](x,z)∇B*P), for ∇⋅j;
* The Boundary Value Problem that is currently implemented assumes Neumann boundary conditions at all boundaries.
Developed by:\\\n
Niels Grobbe, Massachusetts Institute of Technology, USA.\\\n
In collaboration with: Aimé Fournier & Laurent Demanet, Massachusetts Institute of Technology, USA.\\\n
Date: October, 2017 \\\n
Contact: ngrobbe@gmail.com
"""
module Poisson

using SparseArrays
using StatsBase
using LinearAlgebra
using Random
using ProgressMeter
import JuMIT.Data
import JuMIT.Acquisition

const CHUNK_SIZE = 10

mutable struct Paramx{T}
	xin::Vector{T}
	xout::Vector{T}
end

mutable struct Param{T,S}
	nx::Int64
	nz::Int64
	dx::Float64
	dz::Float64
	A::SparseMatrixCSC{Float64,Int64}
	x::Paramx{T}
	dual_x::Paramx{S}
end

function Param(nz,nx)
	dz= inv(nz - 1);
	dx= inv(nx - 1);

	z=(0 : nz - 1)*dz;
	x=(0:nx-1)*dx;
	    
	#--------------------  Build explicit operator matrix A=∇⋅(mpars(x,z)∇) for direct Poisson solver --------------------    
	# Some notes: 
	# We build the operator matrix explicitly like this to properly deal with the Neumann boundary conditions at all boundaries. 
	# The standard Kronecker product approach for building sparse operator matrices seemed to generate some errors at the boundaries.
	# Medium properties calculated on staggered grid.  
	# Directions: z-direction row-wise, x-direction column-wise. 

	# Allocate non-zero operator entries
	k = 5*(nz - 2)*(nx - 2)+ # pentadiagonal operator entry allocation: 5, since using 5-point finite-difference stencil (star-shaped)
	    2*(2*(nz - 2) + 2*(nx - 2)+  # bidiagonal non-corner 4-boundary allocation # for boundaries (2 points only)
	    4); # bidiagonal corner boundary allocation
	# The sparse() function is often a handy way to construct sparse matrices. 
	# It takes as its input a vector I of row indices, a vector J of column indices, and a vector V of nonzero values. 
	# Sparse(I,J,V) constructs a sparse matrix such that S[I[k], J[k]] = V[k].
	iz=Array{Float64,1};
	ix=Array{Float64,1};
	iz=collect(1:nz*nx); # same row: row index=1 --> vectors must have same length [k]
	ix=collect(1:nz*nx); # column indices
	val=zeros(nz*nx); # actual values

	A=sparse(ix,iz,val,nz*nx,nz*nx); # pre-allocate sparse operator matrix A
	xin=zeros(nz*nx)
	xout=zeros(nz*nx)

	pax=Paramx(xin,xout)
	dual_pax=Paramx()

	return Param(nx,nz,dx,dz,A, pax, dual_pax)
end


"""
* `mpars` : medium parameters for forward or inverse mode
"""
function updateA!(pa, mpars)

	ddz = inv(pa.dz*pa.dz); # for simplicity below
	ddx = inv(pa.dx*pa.dx);

	updateAcore!(pa.A,mpars,ddx,ddz,pa.nx, pa.nz)
	boundary!(pa)

end

function updateAcore!(A,mpars,ddx,ddz, nx,nz)
	#-------------- Modify operator A for the interior points (without boundary)---------------------------------------
	for j = 2 : nx - 1 # interior x (column) loop (for amount of internal nodes)
	    k = 2 + nz*(j - 1); # initial interior z for this x. To jump from j (new column), need nz*... to proceed to next column vector k
	    for i = 2 : nz - 1  # interior z (row) loop:
		A[k,k - nz] = .5*(mpars[i,j] + mpars[i,j - 1])*ddx; # moving left (and placing actual staggered 2D grid medium value there)
		A[k,k - 1] = .5*(mpars[i,j] + mpars[i - 1,j])*ddz; # moving up the grid
		A[k,k] = -(mpars[i,j] + .5*(mpars[i,j - 1] + mpars[i,j + 1]))*ddx -
		(mpars[i,j] + .5*(mpars[i - 1,j] + mpars[i + 1,j]))*ddz; # diagonal elements: ∇⋅\medpars(x,z)∇ calculated using 5 mpars points
		A[k,k + 1] = .5*(mpars[i,j] + mpars[i + 1,j])*ddz; # moving down the grid
		A[k,k + nz] = .5*(mpars[i,j]+ mpars[i,j + 1])*ddx; # moving right on grid (need pa.nz values to go to next column)
		k = k + 1; # to next z
	    end # The foregoing A terms are operations by hand ...
	end

end

function solve(field_in,pa,solflag)
	field_out=zeros(size(field_in))
	solve!(field_out,field_in,pa,solflag)
	return field_out
end

function boundary!(pa)

	for j = 2 : pa.nx - 1 # interior x (column) loop:
	    k = 1 + pa.nz*(j - 1);# initial boundary z for this x (k = i + m*(j - 1))
	    for i = [1 pa.nz]  # % north & south boundary z (row) loop:
		o = k - Int(sign(i - .5*pa.nz)); #% sign function determines correct direction to: inward from north or south boundary # to define direction(from boundary down or up)
		pa.A[k,k] = -pa.A[o,k]; #% build in symmetry ... --> trick to get symmetry, the next line is leading.
		pa.A[k,o] = -pa.A[k,k]; #% cancelation implies equation # Neumann: value of psi = same. The operator does not know value. Just knows has to subtract two values --> minus sign
		pa.xin[k] = 0; # 0 to equate boundary values in source term ... (?)
		k = k + pa.nz - 1; #% from north to south boundary
	    end
	end


	# Modify A for the east and west boundaries without corner points
	for j = [1 pa.nx] # west & east boundary x (column) loop:
	    k = 2 + pa.nz*(j - 1); #% initial interior z for this x (k = i + m*(j - 1))
	    for i = 2 : pa.nz - 1 # interior z (row) loop:
		o = k - Int(sign(j - .5*pa.nx))*pa.nz; # inward from west or east boundary
		pa.A[k,k] = -pa.A[o,k]; #% build in symmetry ...
		pa.A[k,o] = -pa.A[k,k]; #% cancelation implies equation
		pa.xin[k]= 0;  #% 0 to equate boundary values
		k = k + 1; #% to next z
	    end
	 end

	# Modify A for the corner points
	for j = [1 pa.nx] #% west & east boundary x (column) loop:
	    k = 1 + pa.nz*(j - 1); #% start on north boundary
	    for i = [1 pa.nz] #% north & south boundary y (row) loop:
	       o = k - Int(sign(i - .5*pa.nz))- #% inward from north or south boundary
		Int(sign(j - .5*pa.nx))*pa.nz; # inward from west or east boundary # two terms to move diagonal
		pa.A[k,k] = pa.A[o,o]; #% nopa.nzero, alas asymmetric ...
		pa.A[k,o] = -pa.A[k,k]; #% ... cancelation implies equation
		pa.xin[k] = 0; #% 0 to equate boundary values 
		k = k + pa.nz - 1; #% from north to south boundary
	    end
	 end
	 return nothing

end

"""
Poisson.solve(field,pa,solflag) solves the forward or inverse Poisson problem in a heterogeneous medium.
# Function input arguments:
* `field` : either source term used in 'forward' mode to get field, or the field (e.g. pore pressure) used in 'inverse' mode to calculate source-term 
* `pa` : Param object, or allocated A
* `solflag` : Solution flag determining whether to use forward (solflag=1) or inverse mode (solflag=-1) of Poisson solver.
"""
function solve!(field_out, field_in, pa,solflag::Int64;A=pa.A) # fields: either field or source


	(abs(solflag)≠1) && error("In Poisson.solve: change solflag to solflag=1 (forward problem) or solflag=-1 (inverse problem)")

	#------------- Modify operator A for the top and bottom boundaries including Neumann conditions (without the corners)-------------

	if(solflag==-1) # for 'inverse' problem: source plays no role here, so just allocate zeros with proper size.
		#src=pa.xin

	elseif(solflag==1) # for 'forward' problem
		#src=pa.xin
		copyto!(pa.xin,field_in)
	end


	for j = 2 : pa.nx - 1 # interior x (column) loop:
	    k = 1 + pa.nz*(j - 1);# initial boundary z for this x (k = i + m*(j - 1))
	    for i = [1 pa.nz]  # % north & south boundary z (row) loop:
		pa.xin[k] = 0; # 0 to equate boundary values in source term ... (?)
		k = k + pa.nz - 1; #% from north to south boundary
	    end
	end


	# Modify A for the east and west boundaries without corner points
	for j = [1 pa.nx] # west & east boundary x (column) loop:
	    k = 2 + pa.nz*(j - 1); #% initial interior z for this x (k = i + m*(j - 1))
	    for i = 2 : pa.nz - 1 # interior z (row) loop:
		pa.xin[k]= 0;  #% 0 to equate boundary values
		k = k + 1; #% to next z
	    end
	 end

	# Modify A for the corner points
	for j = [1 pa.nx] #% west & east boundary x (column) loop:
	    k = 1 + pa.nz*(j - 1); #% start on north boundary
	    for i = [1 pa.nz] #% north & south boundary y (row) loop:
		pa.xin[k] = 0; #% 0 to equate boundary values 
		k = k + pa.nz - 1; #% from north to south boundary
	    end
	 end

	if(solflag == -1) #for inverse problem: fields + medium --> source term
		#pp_vec=pp[:]; # vectorize 2D pore pressure field (according to organization operator A).
		#divj=-reshape(pa.A*pp_vec,pa.nz,pa.nx);

		copyto!(pa.xin, field_in)
		mul!(pa.xout, A, pa.xin)
		copyto!(field_out, pa.xout)
		rmul!(field_out, -1.0)

		#return field_out  # function command --> return these variables, in this case electrical source term ∇⋅j (used in inverse problem)

	elseif(solflag==1) # forward problem: source+medium --> fields
		# To 'fix' the nullspace, and make the source consistent with the nullspace --> not always fast/memory allowable since
		# julia backslash operator does not like non-square matrices... (?)
		#ix=Array{Float64,1};
		#ix=ones(pa.nz*pa.nx); # same row: row index=1
		#iz=collect(1:pa.nz*pa.nx); # column indices
		#v=ones(pa.nz*pa.nx); # actual values
		#onemat=sparse(ix,iz,v);
		# A = [onemat;A]; # fix null space
		# psi = reshape(A\[0; src], pa.nz, pa.nx); # make source consistent with null space


		LinearAlgebra.ldiv!(pa.xout,A,pa.xin) # direct inverse solve using backslash operator
		#pa.xout=pa.A\pa.xin # direct inverse solve using backslash operator
		mpsi=StatsBase.mean(pa.xout)
		for i in eachindex(field_out)
			field_out[i]=pa.xout[i] - mpsi
		end
		#psi = reshape(psi_vec, pa.nz, pa.nx); # reshape electrical potential vector psi into 2D electrical potential map (z,x)
		#psi = psi .- StatsBase.mean(StatsBase.mean(psi)); # alternative to 'fix' nullspace: make zero-mean...
	#	src=reshape(src,pa.nz,pa.nx); # reshape source vector into 2D source map for plotting purposes only.
	    
		#return field_out  # function command --> return these variables, in this case electrical potential ψ
	end # end solflag if statement

	return nothing

end # function command


include("PoissonExpt.jl")

end # end module Poisson
