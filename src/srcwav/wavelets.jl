
"""
Generate a Ricker Wavelet. Reference:
Frequencies of the Ricker wavelet, Yanghua Wang, GEOPHYSICS, VOL. 80, NO. 2
Bandwidth is roughly 1.2 * fpeak for ricker

# Keyword Arguments

* `fqdom::Float64`: dominant frequency 
* `tgrid`: time-domain grid
* `tpeak::Float64=tgrid[1]+1.5/fqdom`: the peak of the ricker in time (has a default)
* `tpeak_rand_flag` : randomly place peak of ricker
"""
function ricker(fqdom::Float64,
		tgrid::StepRangeLen;
		tpeak::Float64=tgrid[1]+1.5/fqdom, # using approximate half width of ricker
		attrib::AbstractString="",
		trim_tol::Float64=0.0,
		maxamp::Float64=1.0,
		tpeak_rand_flag=false,
		maxamp_rand_flag=false,
		phase_rand_flag=false,
		)
	(tpeak < tgrid[1]+1.5/fqdom) && error("cannot output Ricker for given tgrid and tpeak")
	(tpeak > tgrid[end]-1.5/fqdom) && error("cannot output Ricker for given tgrid and tpeak")

	if(tpeak_rand_flag)
		tpeak=rand(Uniform(tgrid[1]+1.5/fqdom,tgrid[end]-1.5/fqdom))
	end
	if(maxamp_rand_flag)
		maxamp=rand(Uniform(-1,1))
	end
	isapprox(fqdom,0.0) && error("dominant frequency cannot be zero")

	#! some constants
	pf = (π*π)*(fqdom^2.0)
	nt = length(tgrid)
	δt = step(tgrid)

	# a vector is odd number of samples (nt + 1 corresponds to time zero)
	wav = zero(tgrid);
	# k = (1 - 2* pf * t^2) * Exp[-pf *t^2]
	# Simplify[D[k,t]]
	# FortranForm[Simplify[D[k,t]]]
	if(occursin("[DIFF]",attrib))
			# ricker after a time derivative
			for it = 1:nt
				tsquare = (tgrid[it]-tpeak) * (tgrid[it]-tpeak)
				t       = -1.0 * (tgrid[it]-tpeak)
				wav[it] = (2.0 * pf * t * (-3.0 + 2.0 * pf * tsquare)) * exp(-1.0 * pf * tsquare)
			end
	else
	#! ricker wavelet
		for it = 1:nt
			tsquare = (tgrid[it]-tpeak) * (tgrid[it]-tpeak)
			wav[it] = (1.0 - 2.0 * pf * tsquare) * exp(-1.0e0 * pf * tsquare) * maxamp
		end
	end

	isapprox(maximum(abs.(wav)),0.0) && warn("wavelet is zeros")

	# apply random phase to the ricker wavelet
	if(phase_rand_flag)
		wav=apply_rand_phase(wav)
	end

	if(trim_tol != 0.0)
		return wav[abs(wav).>=trim_tol]
	else
		return wav
	end
end


"""
ormbsy wavelet


* `tperc::Float64=0.0` : the wavelet is tapered in time using this percentage value
"""
function  ormsby(
		fqdom,
		tgrid::StepRangeLen;
		fracbandwidth::Float64=1.2, # using same value as Ricker
		f1::Float64=fqdom-fracbandwidth*0.5*fqdom,
		f4::Float64=fqdom+fracbandwidth*0.5*fqdom,
		f2::Float64=0.25*f4+0.75*f1,
		f3::Float64=0.25*f1+0.75*f4,
		tpeak::Float64=tgrid[1]+1.5/(0.25*(f1+f4)), # using approximate half width 
		trim_tol::Float64=0.0,
		tperc::Float64=0.0,
		tpeak_rand_flag=false,
		maxamp_rand_flag=false,
		phase_rand_flag=false,
		)

	# some constants
	A43 = (pi*f4)^2 / (pi*f4 - pi*f3);
	A34 = (pi*f3)^2 / (pi*f4 - pi*f3);
	A21 = (pi*f2)^2 / (pi*f2 - pi*f1);
	A12 = (pi*f1)^2 / (pi*f2 - pi*f1);

	wav = zeros(length(tgrid));

	if(tpeak_rand_flag)
		tpeak=rand(Uniform(tgrid[1]+1.5/(0.25*(f1+f4)),tgrid[end]-1.5/(0.25*(f1+f4))))
	end

	# ormsby wavelet
	for it = 1: length(tgrid)
		t = tgrid[it]-tpeak
		S4 = (sinc(pi*f4*t))^2
		S3 = (sinc(pi*f3*t))^2
		S2 = (sinc(pi*f2*t))^2
		S1 = (sinc(pi*f1*t))^2

		wav[it] =  ((A43 * S4 - A34 * S3) - (A21 * S2 - A12 * S1))
	end

	isapprox(maximum(abs.(wav)),0.0) && warn("wavelet is zeros")

	# normalize
	wav /= maximum(abs.(wav));

	if(maxamp_rand_flag)
		wav *= rand(Uniform(-1,1))
	end



	# apply random phase to the ricker wavelet
	if(phase_rand_flag)
		wav=apply_rand_phase(wav)
	end

	if(trim_tol != 0.0)
		wav = wav[abs(wav).>=trim_tol]
	end

	wav = wav .* Utils.taper(ones(length(wav)),tperc)
	return wav

end


function apply_rand_phase(wav)
	nt=length(wav)
	W=FFTW.rfft(wav);
	θ=rand(Uniform(-Float64(pi),Float64(pi)))
	for i in eachindex(W)
		W[i]=W[i]*complex(cos(θ),sin(θ))
	end
	wav=FFTW.irfft(W,nt)
end

"""
Return a ricker wavelet for a given input `Medium`.

# Arguments
* `mod::Medium` : Medium
* `nλ::Int64=10` : number of wavelengths (P-wave) in the medium
* `tmaxfrac::Float64=1.0` : by default the maximum modelling time is computed using the average velocity and the diagonal distance of the medium, 
use this fraction to increase of reduce the maximum time

# Keyword Arguments
* all the keywords arguments of the `ricker` method can be used.
"""
function ricker(mod::Medium, nλ::Int=10, tmaxfrac::Real=1.0, epsilon=inv(sqrt(ndims(mod))); args... )
	@assert(!iszero(mod))
	fqdom, tgrid = get_fqdom_tgrid(mod, nλ, tmaxfrac, epsilon)
	wav=ricker(fqdom, tgrid; args...)
	return wav, tgrid
end

"""
Same as ricker, but return ormsby...
"""
function ormsby(mod::Medium, nλ::Int=10, tmaxfrac::Real=1.0, epsilon=inv(sqrt(ndims(mod))); args... )
	@assert(!iszero(mod))
	fqdom, tgrid = get_fqdom_tgrid(mod, nλ, tmaxfrac, epsilon)
	wav=ormsby(fqdom, tgrid; args...)
	return wav, tgrid
end

"""
Return dominant source frequency, and its temporal grid for a finite-difference simulation, for given number of wavelengths `nλ` in the medium.
The model has `nλ` wavelengths, and the maximum modeling time is determined by `tmaxfrac`.
`epsilon` is Courant number.
"""
function get_fqdom_tgrid( mod::Medium, nλ::Int, tmaxfrac::Real, epsilon)

	# maximum distance (diagnol) the wave travels
	d = sqrt(sum([(m[1]-m[end])^2 for m in mod.mgrid]))

	# dominant wavelength using mod dimensions
	λdom=d*inv(real(nλ))

	# average P velocity
	vavg=mod.ref[:vp]

	fqdom = vavg/λdom

	# use two-way maximum distance to get tmax
	tmax=2.0*d*inv(vavg)*tmaxfrac

	# choose sampling interval to obey max freq of source wavelet
	δmin = minimum(step.(mod.mgrid))
	vmax=try
		sqrt(mod.bounds[:vp][2]^2 + mod.bounds[:vs][2]^2) # see Virieux (1986)
	catch
		mod.bounds[:vp][2]
	end
	δt=epsilon*δmin/vmax

	# check if δt is reasonable
	#(δt > 0.1/fqdom) : error("decrease spatial sampling or nλ")

	tgrid=range(0.0, stop=tmax, step=δt)

	return fqdom, tgrid
end


