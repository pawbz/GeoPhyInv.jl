
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

"""
transpose of dAdx
"""
function dAdmpars(dx,dz,nx,nz,T)

	ddz = inv(dz*dz); # for simplicity below
	ddx = inv(dx*dx);
	X=[spzeros(T,(nz*nx)+1,nz*nx) for i in 1:nz*nx]  # note the size, it has appended zeros

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

		X[imp][k1] += .5*ddx; ; X[imp_0_m1][k1] += .5*ddx; ;
		X[imp][k2] += .5*ddz; ; X[imp_m1_0][k2] += .5*ddz; ;

		X[imp][k3] -= ddx  ; X[imp][k3] -= ddz; ;
		X[imp_0_m1][k3] -= .5*ddx;; X[imp_0_p1][k3] -= .5*ddx;;
		X[imp_m1_0][k3] -= .5*ddz;; X[imp_p1_0][k3] -= .5*ddz;;

		X[imp][k4] += .5*ddz; ; X[imp_p1_0][k4] += .5*ddz; ;
		X[imp][k5] += .5*ddx; ; X[imp_0_p1][k5] += .5*ddx; ;
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
		for chu in 1:nz*nx
			X[chu][k1] = -X[chu][k3]; 
			X[chu][k2] = -X[chu][k1]; 
		end
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
		for chu in 1:nz*nx
			X[chu][k1] = -X[chu][k3]; 
			X[chu][k2] = -X[chu][k1]; 
		end
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
		for chu in 1:nz*nx
			X[chu][k1] = X[chu][k3]; 
			X[chu][k2] = -X[chu][k1]; 
		end
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
	# the average of four corners would go to zero
	A[nz*nx+1,1] = 1.0
	A[nz*nx+1,nz] = 1.0
	A[nz*nx+1,nz*(nx-1)+1] = 1.0
	A[nz*nx+1,nz*nx] = 1.0

	return nothing
end

