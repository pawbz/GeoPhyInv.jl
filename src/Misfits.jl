__precompile__()

module Misfits

import JuMIT.Data 

"""
Normalized least-squares error between two arrays after 
estimating a scalar that best fits on to another.
Return misfit and α such that αx-y is minimum.
Normalization is done with respect to the 
norm of y.
"""
function error_after_scaling{T}(
			     x::Array{T},
			     y::Array{T}
			    )
	(size(x) ≠ size(y)) && error("x and y different sizes") 
	α = sum(x.*y)/sum(x.*x)
	J = norm(y-α*x)/norm(y)

	return J, α
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

	normfact = Data.TD_dot(y, y)
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

function fg_cls_conv(r, s, w)


end


end # module
