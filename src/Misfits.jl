__precompile__()

module Misfits

import JuMIT.Data 

"""
Normalized least-squares error between two vectors after 
estimating a scalar that best fits on to another.
Return misfit and α such that αx-y is minimum.
Normalization is done with respect to the 
norm of y.
"""
function error_after_scaling(
			     x::Array{Float64},
			     y::Array{Float64}
			    )
length(x) != length(y) ? error("vectors of different length") :

α = dot(x,y)/dot(x,x)

J = norm(y-α*x)/norm(y)
#J = norm(α.*x-y)/norm(y)

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
	     w::Data.TD=Data.TD_ones(x.nfield,x.tgrid,x.acqgeom))

	# check if x and y are similar
	tgrid = x.tgrid;
	acq = x.acqgeom;
	nfield = x.nfield;
	nss = acq.nss;

	f = 0.0;
	for ifield=1:x.nfield, iss=1:acq.nss, ir=1:acq.nr[iss]
		if(!(dfdx === nothing))
			ft = fg_cls!(view(dfdx.d[iss,ifield],:,ir), view(x.d[iss, ifield],:,ir), y.d[iss, ifield][:,ir], w.d[iss, ifield][:, ir]);
		else
			ft = fg_cls!(nothing, view(x.d[iss, ifield],:,ir), y.d[iss, ifield][:,ir], w.d[iss, ifield][:, ir]);
		end
		"multiplication with time sampling due to integration ?"
		f += ft #* tgrid.δx;
	end
	f == 0.0 ? warn("misfit is zero") : nothing

	return f
end


function fg_cls!{N}(dfdx, 
		    x::AbstractArray{Float64,N}, 
		    y::Array{Float64,N}, w::Array{Float64,N}=ones(x))
	(size(x) == size(y)) || error("size mismatch")
	any(w .< 0.0) && error("weights cannot be negative") 
	f = sum(w .* (x - y).^2)
	if(!(dfdx === nothing))
		dfdx[:] = 2.0 .* w[:] .* (x[:] - y[:])
	end
	return f
end

function fg_cls_conv(r, s, w)


end


end # module
