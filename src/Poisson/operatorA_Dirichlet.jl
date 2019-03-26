
# useful guide here http://math.mit.edu/~stevenj/18.303/fall13/lecture-10.html

diff1(n, d=1.) = inv(d) .* [ [1.0 zeros(1,n-1)]; diagm(1=>ones(n-1)) - diagm(0=>ones(n)) ]
project_on_diff1(n) = [[1 zeros(1,n-1)]; diagm(1=>0.5.*ones(n-1)) + diagm(0=>vcat(0.5.*ones(n-1), [1.0]))]
project_on_Dx(nz,nx) = kron(project_on_diff1(nx), diagm(0=>ones(nz)))
project_on_Dz(nz,nx) = kron(diagm(0=>ones(nx)), project_on_diff1(nz))
diff1x(nz,nx,dz=1.,dx=1.) = kron(diff1(nx,dx), diagm(0=>ones(nz)))
diff1z(nz,nx,dz=1.,dx=1.) = kron(diagm(0=>ones(nx)), diff1(nz,dz))


#=
x=repeat([1,2,3], 1,4)
xx=project_on_Dz(3,4)*vec(x)



function Laplacian(Nx, Ny)

	dy=inv(Ny-1);
	dx=inv(Nx-1);
   Dx = diff1(Nx) / dx
   Dy = diff1(Ny) / dy
   Ax = Dx' * Dx
   Ay = Dy' * Dy
   return kron(diagm(0=>ones(Ny)), Ax) + kron(Ay, diagm(0=>ones(Nx)))
end


n=5
A=Array(Param(n,n,ones(n*n)).A)[1:n*n,1:n*n]
B=Laplacian(n,n)

x=zeros(n,n)
x[div(n,2)+1,div(n,2)+1]=1.0
x=vec(x)


xx1=reshape(A*x,n,n)
xx2=reshape(B*x,n,n)
=#



"""
* `mpars` : medium parameters for forward or inverse mode
"""
function updateA!(pa::Param, mpars; A=pa.A)
	println("ffffwewe")
	nz=pa.nz; nx=pa.nx; dz=pa.dz; dx=pa.dx
	nznx=nz*nx
	dx=pa.dx
	dz=pa.dz

	Ax=diff1x(nz,nx,dz,dx)' * diagm(0=>project_on_Dx(nz,nx)*vec(mpars))  * diff1x(nz,nx,dz,dx)
	Az=diff1z(nz,nx,dz,dx)' * diagm(0=>project_on_Dz(nz,nx)*vec(mpars))  * diff1z(nz,nx,dz,dx)
	
	copyto!(A,Ax.+Az)

	# update At
	for i in 1:nznx+1
		for j in 1:nznx
			pa.At[j,i]=A[i,j]
		end
	end
	return pa
end


"""
transpose of dAdx
"""
function dAdmpars(dx,dz,nx,nz,T)

	ddz = inv(dz*dz); # for simplicity below
	ddx = inv(dx*dx);
	X=[spzeros(T,(nz*nx)+1,nz*nx) for i in 1:nz*nx]  # note the size, it has appended zeros

	m1=randn(nz*nx)
	m2=zeros(nz*nx)

	for i in 1:nz*nx
		copyto!(m2,m1)
		m2[i]+=0.1
		Ax1=diff1x(nz,nx,dz,dx)' * diagm(0=>project_on_Dx(nz,nx)*m1)  * diff1x(nz,nx,dz,dx)
		Az1=diff1z(nz,nx,dz,dx)' * diagm(0=>project_on_Dz(nz,nx)*m1)  * diff1z(nz,nx,dz,dx)

		Ax2=diff1x(nz,nx,dz,dx)' * diagm(0=>project_on_Dx(nz,nx)*m2)  * diff1x(nz,nx,dz,dx)
		Az2=diff1z(nz,nx,dz,dx)' * diagm(0=>project_on_Dz(nz,nx)*m2)  * diff1z(nz,nx,dz,dx)
		copyto!(X[i],Ax1+Az1-Ax2-Az2)
		rmul!(X[i],10)
	end


	 return X
end


