



"""
Type for inversion variable
* `x` : inversion vector that is optimized
* `m` : model without precon
* `gx` : allocation to store gradient of x 
"""
mutable struct X{T,N}
	x::Vector{T}
	last_x::Vector{T}  # last_x needed for optim
	lower_x::Vector{T}
	upper_x::Vector{T}
	m::Array{T,N}
	gx::Vector{T}
	gm::Vector{Array{T,N}}
	prior::Vector{T} # prior on x
	precon::Vector{T}
	preconI::Vector{T} # inv of precon
	w::Vector{T} # only used to compute the distance from prior
end

function X(nx::Int, nobj::Int, m=nothing)
	x=zeros(nx)
	last_x=zero(x)
	lower_x=zero(x)
	upper_x=zero(x)
	gx=zeros(nx)
	prior=zeros(nx)

	precon=ones(nx)
	preconI=ones(nx)

	w=ones(nx)

	if(m===nothing)
		m=zeros(nx)
	end
	gm=[zero(m) for iobj in 1:nobj]

	xx=X(x,last_x, upper_x, lower_x, m,gx,gm,prior,precon,preconI,w)

	update_preconI!(xx)
	return xx
end



"""
update inversion variable using model
Note that `m` can be an array of length same as `x`.
"""
function update_x!(X)
	for i in eachindex(X.x)
		@inbounds X.x[i]=X.m[i]*X.precon[i] 
	end
end

function update_m!(X)
	for i in eachindex(X.x)
		@inbounds X.m[i]=X.x[i]*X.preconI[i] 
	end
end


"""
update ∂ₓJ using ∂ₘJ
"""
function update_gx!(X, α)
	gx=X.gx
	for i in eachindex(gx)
		gx[i]=zero(eltype(gx))
	end
	for (j,gm) in enumerate(X.gm)
		for i in eachindex(gx)
			@inbounds gx[i]+=gm[i]*X.preconI[i]*α[j]
		end
	end
end

function initialize!(X)
	T=eltype(X.x)
	for i in eachindex(X.x)
		X.x[i]=zero(T)
		X.m[i]=zero(T)
		X.gx[i]=zero(T)
		X.prior[i]=zero(T)
		X.precon[i]=one(T)
		X.preconI[i]=one(T)
		X.w[i]=one(T)
	end
	for j in 1:length(X.gm)
		X.gm[j][:]=zero(T)
	end
end


function update_precon!(X, precon)
	p=X.precon
	if(precon===nothing)
		for i in eachindex(p)
			p[i]=one(p)
		end
	else
		copyto!(p, precon)
	end
	update_preconI(X)
end
function update_preconI!(X)
	for i in eachindex(X.precon)
		if(iszero(X.precon[i]))
			X.preconI[i]=X.precon[i]
		else
			X.preconI[i]=inv(X.precon[i])
		end
	end
end

function update_weights!(X, weights)
	w=X.w
	if(weights===nothing)
		for i in eachindex(w)
			w[i]=one(eltype(w))
		end
	else
		copyto!(w, weights)
	end
end


