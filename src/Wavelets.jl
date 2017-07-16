__precompile__()

module Wavelets

import SIT.Grid

"""
Generate a Ricker Wavelet. Reference:
Frequencies of the Ricker wavelet, Yanghua Wang, GEOPHYSICS, VOL. 80, NO. 2

# Keyword Arguments

* `fqdom::Float64`: dominant frequency 
* `tgrid::Grid.M1D`: time-domain grid
* `tpeak::Float64=tgrid.x[1]+1.5/fqdom`: the peak of the ricker in time (has a default)
"""
function ricker(;
		fqdom::Float64=nothing,
		tgrid::Grid.M1D=nothing,
		tpeak::Float64=tgrid.x[1]+1.5/fqdom, # using approximate half width of ricker
		attrib::AbstractString="",
		trim_tol::Float64=0.0
		)
	(tpeak < tgrid.x[1]+1.5/fqdom) ? error("cannot output Ricker for given tgrid and tpeak") : nothing
	(tpeak > tgrid.x[end]-1.5/fqdom) ? error("cannot output Ricker for given tgrid and tpeak")  : nothing

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

	isapprox(maximum(abs.(wav)),0.0) && warn("wavelet is zeros")

	if(trim_tol != 0.0)
		return wav[abs(wav).>=trim_tol]
	else
		return wav
	end
end


"""
ormbsy wavelet
"""
function  ormsby(;
		f1::Float64=5.0,
		f4::Float64=40.0,
		f2::Float64=0.25*f4+0.75*f1,
		f3::Float64=0.25*f1+0.75*f4,
		tpeak::Float64=0.25,
		tgrid::Grid.M1D=nothing,
		trim_tol::Float64=0.0
		)

# some constants
A43 = (pi*f4)^2 / (pi*f4 - pi*f3);
A34 = (pi*f3)^2 / (pi*f4 - pi*f3);
A21 = (pi*f2)^2 / (pi*f2 - pi*f1);
A12 = (pi*f1)^2 / (pi*f2 - pi*f1);

wav = zeros(tgrid.nx);
# ormsby wavelet
for it = 1: tgrid.nx
	t = tgrid.x[it]-tpeak
        S4 = (sinc(pi*f4*t))^2
        S3 = (sinc(pi*f3*t))^2
        S2 = (sinc(pi*f2*t))^2
        S1 = (sinc(pi*f1*t))^2

	wav[it] =  (A43 * S4 - A34 * S3) - (A21 * S2 - A12 * S1)
end

isapprox(maximum(abs(wav)),0.0) && warn("wavelet is zeros")

# normalize
wav /= maximum(abs(wav));

if(trim_tol != 0.0)
	return wav[abs(wav).>=trim_tol]
else
	return wav
end

end

end # module
