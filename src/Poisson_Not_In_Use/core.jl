using SparseArrays
using StatsBase
using LinearAlgebra
using Random
using ProgressMeter
using ForwardDiff


mutable struct Param{T}
	A::SparseMatrixCSC{T,Int64}
	At::SparseMatrixCSC{T,Int64}  # for adjoint solver
	dAdx::Vector{SparseMatrixCSC{T,Int64}}  # necessary for the gradient
	x1::Vector{T}
	x2::Vector{T}
	nx::Int64
	nz::Int64
	dx::Float64
	dz::Float64
end

include("operatorA_Neumann.jl")
#include("operatorA_Dirichlet.jl")

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

	dAdx=dAdmpars(dx,dz,nx,nz,T)

	pa=Param(A,At,dAdx,zeros(T, nznx),zeros(T, nznx+1),nx,nz,dx,dz)

	if(!(mpars===nothing))
		updateA!(pa, mpars)
	end

	return pa
end



"""
"""
function applyinvA!(fout, fin, pa; A=pa.A, mute_boundary_source=true) 
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

	if(mute_boundary_source)
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
	fill!(fout, 0.0) # clear output
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
