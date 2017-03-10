module Wavelets

import SIT.Grid: M1D

function ricker(;
		fqdom::Float64=20.0,
		tgrid::M1D=M1D(:time),
		tpeak::Float64=0.25,
		attrib::AbstractString="",
		trim_tol::Float64=0.0
		)

isapprox(fqdom,0.0) && error("dominant frequency cannot be zero")

#! some constants
pf = (π*π)*(fqdom^2.0)
nt = tgrid.nx
δt = tgrid.δx

# a vector is odd number of samples (nt + 1 corresponds to time zero)
wav = zeros(tgrid.x);
# k = (1 - 2* pf * t^2) * Exp[-pf *t^2]
# Simplify[D[k,t]]
# FortranForm[Simplify[D[k,t]]]
if(contains(attrib,"[DIFF]"))
                # ricker after a time derivative
                for it = 1:nt
			tsquare = (tgrid.x[it]-tpeak) * (tgrid.x[it]-tpeak)
			t       = -1.0 * (tgrid.x[it]-tpeak)
                        wav[it] = (2.0 * pf * t * (-3.0 + 2.0 * pf * tsquare)) * exp(-1.0 * pf * tsquare)
                end
else
#! ricker wavelet
        for it = 1:nt
		tsquare = (tgrid.x[it]-tpeak) * (tgrid.x[it]-tpeak)
                wav[it] = (1.0 - 2.0 * pf * tsquare) * exp(-1.0e0 * pf * tsquare)
        end
end

isapprox(maximum(abs(wav)),0.0) && warn("wavelet is zeros")

if(trim_tol != 0.0)
	return wav[abs(wav).>=trim_tol]
else
	return wav
end
end

end # module
