__precompile__()

module Misfits

using ForwardDiff
using Distances
using StatsBase
import JuMIT.Data 
import JuMIT.Grid


"""
Distance between x[:,i] and x[:,j]
cosine_dist(x, y)=1 - dot(x, y) / (norm(x) * norm(y))	
corr_dist(x, y)=cosine_dist(x - mean(x), y - mean(y))	
"""
function error_pairwise_corr_dist(dfdx,x)
	R=pairwise(CorrDist(), x)
	J=sum(R)

	if(!(dfdx===nothing))
		f(x)=error_pairwise_corr_dist(nothing,x)
		ForwardDiff.gradient!(dfdx, f, x);
	end
	return J
end

function error_corr_dist!(dfdx, x, y)
	xn=vecnorm(x)
	yn=vecnorm(y)

	dotxy=0.0
	for i in eachindex(x)
		dotxy+=x[i]*y[i]
	end

	J = 1. - (dotxy/xn/yn)
	
	if(!(dfdx===nothing))
		for i in eachindex(x)
			dfdx[i]=(dotxy*yn*x[i]/xn - y[i]*yn*xn)/(xn*yn)/(xn*yn)
		end
	end
	return J
end


function error_autocorr_pairwise_corr_dist(dfdx,x)
	X=fft(x, [1]);
	# remove phase
	#XC=real(X .* conj(X))
	J=error_pairwise_corr_dist(nothing, X)
	J=norm(J)

	if(!(dfdx===nothing))
		GXC=similar(X)
		error_pairwise_corr_dist(GXC,X)
		#GX = GXC .* conj(X)
		dfdx[:]=real.(ifft(GXC, [1]))
	end

	return J
end


function error_weighted_norm(dfdx,x,w)
	J=0.0
	for i in eachindex(x)
		J += (x[i]*w[i])*(x[i]*w[i])
	end
	J=sqrt(J)
	if(!(dfdx===nothing))
		f(x)=error_weighted_norm(nothing,x,w)
		ForwardDiff.gradient!(dfdx, f, x);
	end
	return J
end



"""
Compute the error with a normalized vector ̂̂x 
̂x=x/|x|

function 
* `ghat` : gradient with respect to ̂x
* `g`	 : output gradient with respect to x
* `X`	 : preallocated matrix of dimension (nx,nx) if necessary
"""
function derivative_vector_magnitude!(g,ghat,x,X=nothing) 
	xn=vecnorm(x)
	nx=length(x)
	scale!(x, inv(xn))  

	# compute the outer product of 
	if(!(X===nothing))      
		A_mul_Bt!(X,x,x)        
	else    
		X=A_mul_Bt(x,x)
	end         
	for i in 1:nx
		X[i,i]=X[i,i]-1. 
	end
	scale!(X,-inv(xn))
	A_mul_B!(g,X,ghat)    
	scale!(x, xn) 
end

"""
Normalized least-squares error between two arrays after 
estimating a scalar that best fits on to another.
Return misfit and α such that αx-y is minimum.
Normalization is done with respect to the 
norm of y.
"""
function error_after_scaling{T}(
			     x::AbstractArray{T},
			     y::AbstractArray{T}
			    )
	any(size(x) ≠ size(y)) && error("x and y different sizes") 
	α = sum(x.*y)/sum(x.*x)
	J = norm(y-α*x)/norm(y)

	return J, α
end


"""
Measure the least-squares distance between auto correlations of x and y, 
"""
function error_after_normalized_autocor(x::AbstractArray{Float64}, y::AbstractArray{Float64})
	nt=size(x,1)

	ax=autocor(x, 0:nt-1, demean=true) # normalizd autocorr from StatsBase
	ay=autocor(y, 0:nt-1, demean=true)

	return error_squared_euclidean!(nothing,  ax,   ay, nothing, norm_flag=true)
end

"""
Input the obeserved and modelled data to output the misfit
and the adjoint sources
TODO: 

"""
function TD!(dfdx,
	     x::Data.TD, 
	     y::Data.TD, 
	     w::Data.TD=Data.TD_ones(x.fields,x.tgrid,x.acqgeom))

	# check if x and y are similar
	tgrid = x.tgrid;
	acq = x.acqgeom;
	fields = x.fields;
	nss = acq.nss;

	normfact = dot(y, y)
	(normfact == 0.0) && error("y cannot be zero")

	f = 0.0;
	for ifield=1:length(fields), iss=1:acq.nss, ir=1:acq.nr[iss]
		if(!(dfdx === nothing))
			ft = fg_cls!(view(dfdx.d[iss,ifield],:,ir), view(x.d[iss, ifield],:,ir), y.d[iss, ifield][:,ir], w.d[iss, ifield][:, ir]);
			dfdx.d[iss,ifield][:,ir] /= normfact
		else
			ft = fg_cls!(nothing, view(x.d[iss, ifield],:,ir), y.d[iss, ifield][:,ir], w.d[iss, ifield][:, ir]);
		end
		"multiplication with time sampling due to integration ?"
		f += ft / normfact #* tgrid.δx;
	end
	(f == 0.0) && warn("misfit computed is zero")

	return f
end


function fg_cls!{N}(dfdx, 
		    x::AbstractArray{Float64,N}, 
		    y::Array{Float64,N}, 
		    w::Array{Float64,N}=ones(x))
	(size(x) == size(y) == size(w)) || error("sizes mismatch")
	f = sum(w .* (x - y).^2)
	if(!(dfdx === nothing))
		copy!(dfdx, 2.0 .* w .* (x-y))
	end
	return f
end


function error_squared_euclidean!(dfdx,  x,   y,   w; norm_flag=false)
	J=zero(Float64)
	if(w===nothing)
		for i in eachindex(x)
			J += (x[i]-y[i]) * (x[i]-y[i])
		end
	else
		for i in eachindex(x)
			J += w[i] * (x[i]-y[i]) * (x[i]-y[i])
		end
	end
	if(!(dfdx === nothing))
		if(w===nothing)
			for i in eachindex(x)
				dfdx[i] = 2.0 * (x[i]-y[i])
			end
		else
			for i in eachindex(x)
				dfdx[i] = 2.0 * w[i] * (x[i]-y[i])
			end
		end
	end
	if(norm_flag) # second norm of y
		ynorm=zero(Float64)
		for i in eachindex(y)
			ynorm += y[i]*y[i]
		end
		# divide the functional with ynorm
		J /= ynorm
		if(!(dfdx === nothing)) 
			scale!(dfdx, inv(ynorm)) # divide dfdx with ynorm too
		end
	end
	return J
end


function error_weighted_norm!(dfdx,  err,   w)
	J=zero(Float64)
	if(w===nothing)
		for i in eachindex(err)
			J += (err[i]) * (err[i])
		end
	else
		for i in eachindex(err)
			J += w[i] * (err[i]) * (err[i])
		end
	end
	if(!(dfdx === nothing))
		if(w===nothing)
			for i in eachindex(err)
				dfdx[i] = 2.0 * (err[i])
			end
		else
			for i in eachindex(err)
				dfdx[i] = 2.0 * w[i] * (err[i])
			end
		end
	end
	return J
end



function fg_cls_conv(r, s, w)


end


end # module
