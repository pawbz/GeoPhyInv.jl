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
		diff = x.d[iss, ifield][:,ir] - y.d[iss, ifield][:,ir];
		f += sum((diff).^2)
		δx.d[iss, ifield][:,ir] = 2.*diff
	end

	# check for zeros
	Data.TD_iszero(δx) ? error("δx cannot be zero") : nothing

	return f, δx
end


end # module
