__precompile__()

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


"""
Poisson.solve(field,mpars,solflag) solves the forward or inverse Poisson problem in a heterogeneous medium.
# Function input arguments:
* `field` : either source term used in 'forward' mode to get field, or the field (e.g. pore pressure) used in 'inverse' mode to calculate source-term 
* `mpars` : medium parameters for forward or inverse mode
* `solflag` : Solution flag determining whether to use forward (solflag=1) or inverse mode (solflag=-1) of Poisson solver.
"""

function solve(field,mpars,solflag::Int64) # fields: either field or source, mpars=medium pars

if(solflag==-1) 
	pp=field; #
elseif(solflag==1)
	src=field;
else
	error("In Poisson.solve: change solflag to solflag=1 (forward problem) or solflag=-1 (inverse problem)")
end


nz=size(mpars,1);
nx=size(mpars,2);
dz= 1./(nz - 1);
dx= 1./(nx - 1);

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

dz = dz.^2; # for simplicity below
dx = dx.^2;

#-------------- Create operator A for the interior points (without boundary)---------------------------------------
for j = 2 : nx - 1 # interior x (column) loop (for amount of internal nodes)
    k = 2 + nz*(j - 1); # initial interior z for this x. To jump from j (new column), need nz*... to proceed to next column vector k
    for i = 2 : nz - 1  # interior z (row) loop:
        A[k,k - nz] = .5*(mpars[i,j] + mpars[i,j - 1])/dx; # moving left (and placing actual staggered 2D grid medium value there)
        A[k,k - 1] = .5*(mpars[i,j] + mpars[i - 1,j])/dz; # moving up the grid
        A[k,k] = -(mpars[i,j] + .5*(mpars[i,j - 1] + mpars[i,j + 1]))/dx -
        (mpars[i,j] + .5*(mpars[i - 1,j] + mpars[i + 1,j]))/dz; # diagonal elements: ∇⋅\medpars(x,z)∇ calculated using 5 mpars points
        A[k,k + 1] = .5*(mpars[i,j] + mpars[i + 1,j])/dz; # moving down the grid
        A[k,k + nz] = .5*(mpars[i,j]+ mpars[i,j + 1])/dx; # moving right on grid (need nz values to go to next column)
        k = k + 1; # to next z
    end # The foregoing A terms are operations by hand ...
end

#------------- Create operator A for the top and bottom boundaries including Neumann conditions (without the corners)-------------

if(solflag==-1) # for 'inverse' problem: source plays no role here, so just allocate zeros with proper size.
	src=zeros(nz,nx); 

elseif(solflag==1) # for 'forward' problem
	src = src[:]; # vectorize source (which was function input argument)
end



 for j = 2 : nx - 1 # interior x (column) loop:
    k = 1 + nz*(j - 1);# initial boundary z for this x (k = i + m*(j - 1))
    for i = [1 nz]  # % north & south boundary z (row) loop:
        o = k - Int(sign(i - .5*nz)); #% sign function determines correct direction to: inward from north or south boundary # to define direction(from boundary down or up)
        A[k,k] = -A[o,k]; #% build in symmetry ... --> trick to get symmetry, the next line is leading.
        A[k,o] = -A[k,k]; #% cancelation implies equation # Neumann: value of psi = same. The operator does not know value. Just knows has to subtract two values --> minus sign
        src[k] = 0; # 0 to equate boundary values in source term ... (?)
        k = k + nz - 1; #% from north to south boundary
    end
 end


# Create A for the east and west boundaries without corner points
for j = [1 nx] # west & east boundary x (column) loop:
    k = 2 + nz*(j - 1); #% initial interior z for this x (k = i + m*(j - 1))
    for i = 2 : nz - 1 # interior z (row) loop:
        o = k - Int(sign(j - .5*nx))*nz; # inward from west or east boundary
        A[k,k] = -A[o,k]; #% build in symmetry ...
        A[k,o] = -A[k,k]; #% cancelation implies equation
        src[k]= 0;  #% 0 to equate boundary values
        k = k + 1; #% to next z
    end
 end

# Create A for the corner points
for j = [1 nx] #% west & east boundary x (column) loop:
    k = 1 + nz*(j - 1); #% start on north boundary
    for i = [1 nz] #% north & south boundary y (row) loop:
       o = k - Int(sign(i - .5*nz))- #% inward from north or south boundary
        Int(sign(j - .5*nx))*nz; # inward from west or east boundary # two terms to move diagonal
        A[k,k] = A[o,o]; #% nonzero, alas asymmetric ...
        A[k,o] = -A[k,k]; #% ... cancelation implies equation
        src[k] = 0; #% 0 to equate boundary values 
        k = k + nz - 1; #% from north to south boundary
    end
 end

if(solflag == -1) #for inverse problem: fields + medium --> source term
	pp_vec=pp[:]; # vectorize 2D pore pressure field (according to organization operator A).
	divj=-reshape(A*pp_vec,nz,nx);

    	return divj  # function command --> return these variables, in this case electrical source term ∇⋅j (used in inverse problem)



elseif(solflag==1) # forward problem: source+medium --> fields
# To 'fix' the nullspace, and make the source consistent with the nullspace --> not always fast/memory allowable since
# julia backslash operator does not like non-square matrices... (?)
	#ix=Array{Float64,1};
	#ix=ones(nz*nx); # same row: row index=1
	#iz=collect(1:nz*nx); # column indices
	#v=ones(nz*nx); # actual values
	#onemat=sparse(ix,iz,v);
	# A = [onemat;A]; # fix null space
	# psi = reshape(A\[0; src], nz, nx); # make source consistent with null space


	psi_vec=A\src[:]; # direct inverse solve using backslash operator
	psi = reshape(psi_vec, nz, nx); # reshape electrical potential vector psi into 2D electrical potential map (z,x)
	psi = psi - mean(mean(psi)); # alternative to 'fix' nullspace: make zero-mean...
#	src=reshape(src,nz,nx); # reshape source vector into 2D source map for plotting purposes only.
    
    	return psi  # function command --> return these variables, in this case electrical potential ψ
end # end solflag if statement

end # function command



end # end module Poisson
