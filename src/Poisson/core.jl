using SparseArrays
using StatsBase
using LinearAlgebra
using Random
using ProgressMeter
using ForwardDiff

mutable struct Param{T}
	A::SparseMatrixCSC{T,Int64}
	At::SparseMatrixCSC{T,Int64}  # for adjoint solver
	dAdx::SparseMatrixCSC{T,Int64}  # necessary for the gradient
	x1::Vector{T}
	x2::Vector{T}
	nx::Int64
	nz::Int64
	dx::Float64
	dz::Float64
end

function Param(nx::Int,nz::Int,mpars=nothing)

	dz=inv(nz-1);
	dx=inv(nx-1);
	z=(0:nz-1)*dz;
	x=(0:nx-1)*dx;

	return Param([z,x], mpars)
end


function Param(mgrid, mpars=nothing)
	dz=Float64(step(mgrid[1]))
	dx=Float64(step(mgrid[2]))
	nz,nx=length.(mgrid)
		nznx=nz*nx

	if(!(mpars===nothing))
		T=eltype(mpars)
	else
		T=Float64
	end
	A=spzeros(T, nznx+1,nznx)
	At=spzeros(T, nznx,nznx+1)
	dAdx=spzeros(T,nz*nx*(nz*nx+1),nz*nx)  # note the size, it has appended zeros

	pa=Param(A,At,dAdx,zeros(T, nznx),zeros(T, nznx+1),nx,nz,dx,dz)

	if(!(mpars===nothing))
		updateA!(pa, mpars)
	end

	dAdmpars!(pa)
	return pa
end

"""
* `mpars` : medium parameters for forward or inverse mode
"""
function updateA!(pa::Param, mpars; A=pa.A)
	nz=pa.nz; nx=pa.nx; dz=pa.dz; dx=pa.dx
	nznx=nz*nx
	ddz = inv(pa.dz*pa.dz); # for simplicity below
	ddx = inv(pa.dx*pa.dx);
	updateAcore!(A,mpars,ddx,ddz,nx,nz)
	boundary!(A,nx,nz)

	# update At
	for i in 1:nznx+1
		for j in 1:nznx
			pa.At[j,i]=A[i,j]
		end
	end
	return pa
end

"""
#--------------------  Build explicit operator matrix A=∇⋅(mpars(x,z)∇) for direct Poisson solver --------------------    
# Some notes: 
# We build the operator matrix explicitly like this to properly deal with the Neumann boundary conditions at all boundaries. 
# The standard Kronecker product approach for building sparse operator matrices seemed to generate some errors at the boundaries.
# Medium properties calculated on staggered grid.  
# Directions: z-direction row-wise, x-direction column-wise. 

#=
# Allocate non-zero operator entries
k = 5*(nz - 2)*(nx - 2)+ # pentadiagonal operator entry allocation: 5, since using 5-point finite-difference stencil (star-shaped)
    2*(2*(nz - 2) + 2*(nx - 2)+  # bidiagonal non-corner 4-boundary allocation # for boundaries (2 points only)
    4); # bidiagonal corner boundary allocation
# The sparse() function is often a handy way to construct sparse matrices. 
# It takes as its input a vector I of row indices, a vector J of column indices, and a vector V of nonzero values. 

avoid using this subroutine, use `updateA!` instead
"""
function updateAcore!(A::AbstractArray,mpars::AbstractArray,ddx,ddz,nx,nz)
	(length(mpars) ≠ (nz*nx))  && error("error in mpars")
	#-------------- Modify operator A for the interior points (without boundary)---------------------------------------
	for j = 2 : nx - 1 # interior x (column) loop (for amount of internal nodes)
	    k = 2 + nz*(j - 1); # initial interior z for this x. To jump from j (new column), need nz*... to proceed to next column vector k
	    for i = 2 : nz - 1  # interior z (row) loop:
		A[k,k - nz] = .5*(mpars[i+(j-1)*nz] + mpars[i+(j-2)*nz])*ddx; # moving left (and placing actual staggered 2D grid medium value there)
		A[k,k - 1] = .5*(mpars[i+(j-1)*nz] + mpars[i-1+(j-1)*nz])*ddz; # moving up the grid
		A[k,k] = -(mpars[i+(j-1)*nz] + .5*(mpars[i+(j-2)*nz] + mpars[i+(j)*nz]))*ddx -
		(mpars[i+(j-1)*nz] + .5*(mpars[i-1+(j-1)*nz] + mpars[i+1+(j-1)*nz]))*ddz; # diagonal elements: ∇⋅\medpars(x,z)∇ calculated using 5 mpars points
		A[k,k + 1] = .5*(mpars[i+(j-1)*nz] + mpars[i+1+(j-1)*nz])*ddz; # moving down the grid
		A[k,k + nz] = .5*(mpars[i+(j-1)*nz]+ mpars[i+(j)*nz])*ddx; # moving right on grid (need nz values to go to next column)
		k = k + 1; # to next z
	    end # The foregoing A terms are operations by hand ...
	end
	return A

end

function dAdmpars!(pa::Param)
	nz=pa.nz; nx=pa.nx; dz=pa.dz; dx=pa.dx
	ddz = inv(pa.dz*pa.dz); # for simplicity below
	ddx = inv(pa.dx*pa.dx);
	return dAdmpars!(pa.dAdx,ddx,ddz,nx,nz)
end

function dAdmpars!(X,ddx,ddz,nx,nz)

	nznx1=(nz*nx)+1 

	for j = 2 : nx - 1 
	    k = 2 + nz*(j - 1); 
	    for i = 2 : nz - 1  

		imp=i+(j-1)*nz; imp_0_m1=i+(j-2)*nz;
		imp_m1_0=i-1+(j-1)*nz; imp_0_p1=i+(j)*nz;
		imp_p1_0=i+1+(j-1)*nz;

		k1=k+(k-nz-1)*nznx1; k2=k+(k-1-1)*nznx1;
		k3=k+(k-1)*nznx1; k4=k+(k)*nznx1;
		k5=k+(k+nz-1)*nznx1;

		X[k1,imp] += .5*ddx; ; X[k1,imp_0_m1] += .5*ddx; ;
		X[k2,imp] += .5*ddz; ; X[k2,imp_m1_0] += .5*ddz; ;

		X[k3,imp] -= ddx  ; X[k3,imp] -= ddz; ;
		X[k3,imp_0_m1] -= .5*ddx;; X[k3,imp_0_p1] -= .5*ddx;;
		X[k3,imp_m1_0] -= .5*ddz;; X[k3,imp_p1_0] -= .5*ddz;;

		X[k4,imp] += .5*ddz; ; X[k4,imp_p1_0] += .5*ddz; ;
		X[k5,imp] += .5*ddx; ; X[k5,imp_0_p1] += .5*ddx; ;
		k = k + 1; 
	    end 
	end


	for j = 2 : nx - 1 
	    k = 1 + nz*(j - 1);
	    for i = [1 nz]  
		o = k - Int(sign(i - .5*nz)); 
		k1=k+(k-1)*nznx1
		k2=k+(o-1)*nznx1
		k3=o+(k-1)*nznx1
		X[k1,:] .= -X[k3,:]; 
		X[k2,:] .= -X[k1,:]; 
		k = k + nz - 1; 
	    end
	end

	for j = [1 nx] 
	    k = 2 + nz*(j - 1); 
	    for i = 2 : nz - 1 
		o = k - Int(sign(j - .5*nx))*nz; 
		k1=k+(k-1)*nznx1
		k2=k+(o-1)*nznx1
		k3=o+(k-1)*nznx1
		X[k1,:] .= -X[k3,:]; 
		X[k2,:] .= -X[k1,:]; 
		k = k + 1; 
	    end
	 end

	for j = [1 nx] 
	    k = 1 + nz*(j - 1); 
	    for i = [1 nz] 
	       o = k - Int(sign(i - .5*nz))- 
		Int(sign(j - .5*nx))*nz; 
		k1=k+(k-1)*nznx1
		k2=k+(o-1)*nznx1
		k3=o+(o-1)*nznx1
		X[k1,:] .= X[k3,:]; 
		X[k2,:] .= -X[k1,:]; 
		k = k + nz - 1; #
	    end
	 end

	return X
end

function boundary!(A, nx, nz)

	#------------- Modify operator A for the top and bottom boundaries including Neumann conditions (without the corners)-------------
	for j = 2 : nx - 1 # interior x (column) loop:
	    k = 1 + nz*(j - 1);# initial boundary z for this x (k = i + m*(j - 1))
	    for i = [1 nz]  # % north & south boundary z (row) loop:
		o = k - Int(sign(i - .5*nz)); #% sign function determines correct direction to: inward from north or south boundary # to define direction(from boundary down or up)
		A[k,k] = -A[o,k]; #% build in symmetry ... --> trick to get symmetry, the next line is leading.
		A[k,o] = -A[k,k]; #% cancelation implies equation # Neumann: value of psi = same. The operator does not know value. Just knows has to subtract two values --> minus sign
		k = k + nz - 1; #% from north to south boundary
	    end
	end

	# Modify A for the east and west boundaries without corner points
	for j = [1 nx] # west & east boundary x (column) loop:
	    k = 2 + nz*(j - 1); #% initial interior z for this x (k = i + m*(j - 1))
	    for i = 2 : nz - 1 # interior z (row) loop:
		o = k - Int(sign(j - .5*nx))*nz; # inward from west or east boundary
		A[k,k] = -A[o,k]; #% build in symmetry ...
		A[k,o] = -A[k,k]; #% cancelation implies equation
		k = k + 1; #% to next z
	    end
	 end

	# Modify A for the corner points
	for j = [1 nx] #% west & east boundary x (column) loop:
	    k = 1 + nz*(j - 1); #% start on north boundary
	    for i = [1 nz] #% north & south boundary y (row) loop:
	       o = k - Int(sign(i - .5*nz))- #% inward from north or south boundary
		Int(sign(j - .5*nx))*nz; # inward from west or east boundary # two terms to move diagonal
		A[k,k] = A[o,o]; #% nonzero, alas asymmetric ...
		A[k,o] = -A[k,k]; #% ... cancelation implies equation
		k = k + nz - 1; #% from north to south boundary
	    end
	 end

	# append A with ones, to remove mean 
	# To 'fix' the nullspace, and make the source consistent with the nullspace
	#A[nz*nx+1,:] .= 1.0
	A[nz*nx+1,:] .= 0.0
	A[nz*nx+1,1] = 1.0

	return nothing
end

"""
"""
function applyinvA!(fout, fin, pa; A=pa.A) 
	nx=pa.nx
	nz=pa.nz
	nznx=nz*nx

	x1=pa.x1
	x2=pa.x2
	fill!(pa.x1, 0.0)
	fill!(pa.x2, 0.0)

	for i in 1:nznx
		x2[i]=fin[i]
	end

	for j = 2 : nx - 1 # interior x (column) loop:
	    k = 1 + nz*(j - 1);# initial boundary z for this x (k = i + m*(j - 1))
	    for i = [1 nz]  # % north & south boundary z (row) loop:
		x2[k] = 0; # 0 to equate boundary values in source term ... (?)
		k = k + nz - 1; #% from north to south boundary
	    end
	end


	for j = [1 nx] # west & east boundary x (column) loop:
	    k = 2 + nz*(j - 1); #% initial interior z for this x (k = i + m*(j - 1))
	    for i = 2 : nz - 1 # interior z (row) loop:
		x2[k]= 0;  #% 0 to equate boundary values
		k = k + 1; #% to next z
	    end
	 end

	for j = [1 nx] #% west & east boundary x (column) loop:
	    k = 1 + nz*(j - 1); #% start on north boundary
	    for i = [1 nz] #% north & south boundary y (row) loop:
		x2[k] = 0; #% 0 to equate boundary values 
		k = k + nz - 1; #% from north to south boundary
	    end
	 end

#	LinearAlgebra.ldiv!(x1,A,x2) # direct inverse solve using backslash operator
	x1 .= A\x2 # memory allocated here?
	#x11 = A\x2 # memory allocated here?

	for i in 1:nznx
		fout[i]=x1[i]
	end

	return fout
end 

function applyinvAt!(fout, fin, pa; At=pa.At) 
	nx=pa.nx
	nz=pa.nz
	nznx=nz*nx

	x1=pa.x1
	x2=pa.x2
	fill!(pa.x1, 0.0)
	fill!(pa.x2, 0.0)

	for i in 1:nznx
		x1[i]=fin[i]
	end

#	LinearAlgebra.ldiv!(x1,A,x2) # direct inverse solve using backslash operator
	x2 .= At\x1 # memory allocated here?

	for i in 1:nznx+1
		fout[i]=x2[i]
	end

	return fout
end 


"""
Apply `A` to `fin` without allocating memory.
"""
function applyA!(fout, fin, pa; A=pa.A) 
	nx=pa.nx
	nz=pa.nz
	nznx=nz*nx

	x1=pa.x1
	x2=pa.x2
	fill!(pa.x1, 0.0)
	fill!(pa.x2, 0.0)
	for i in 1:nznx
		x1[i]=fin[i]
	end

	mul!(x2, A, x1)
	for i in 1:nznx
		fout[i]=x2[i]
	end
	rmul!(fout, -1.0)
	return fout
end 

function applyinvA(fin,pa;A=pa.A)
	fout=zeros(size(A,2))
	applyinvA!(fout,fin,pa,A=A)
	return fout
end
function applyinvAt(fin,pa;At=pa.At)
	fout=zeros(size(At,2))
	applyinvAt!(fout,fin,pa,At=At)
	return fout
end


function applyA(fin,pa;A=pa.A)
	fout=zeros(size(fin))
	applyA!(fout,fin,pa,A=A)
	return fout
end


