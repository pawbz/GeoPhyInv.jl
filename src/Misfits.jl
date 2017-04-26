module Misfits

import SIT.Data 

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
"""
function TD(x::Data.TD, y::Data.TD
	   )

	# check if x and y are similar
	tgrid = x.tgrid;
	acq = x.acqgeom;
	nfield = x.nfield;
	nss = acq.nss;

	δx=Data.TD([zeros(tgrid.nx,acq.nr[iss]) for iss=1:nss, ifield=1:nfield],
	      nfield,tgrid,acq)

	f = 0.0;
	for ifield=1:x.nfield, iss=1:acq.nss, ir=1:acq.nr[iss]
		ft, δx.d[iss, ifield][:,ir] = fg_cls(x.d[iss, ifield][:,ir], y.d[iss, ifield][:,ir]);
		f += ft;
	end
	f == 0.0 ? warn("misfit is zero") : nothing

	# check for zeros
	Data.TD_iszero(δx) ? error("δx cannot be zero") : nothing

	return f, δx
end


function fg_cls{N}(x::Array{Float64,N}, y::Array{Float64,N})
	size(x) == size(y) ? nothing : error("size mismatch")
	diff = x - y
	f = sum((vec(diff)).^2)
	δx = 2.0 .* diff
	return f, δx
end

function fg_cls_conv(r, s, w)


end


end # module
