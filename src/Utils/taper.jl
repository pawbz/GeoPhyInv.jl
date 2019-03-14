

"""
Cosine taper a N-dimensional array along its first dimension.

# Arguments
* `x::Array{Float64,N}` : 
* `perc::Float64` : taper percentage
"""
function taper(x::AbstractArray,perc::Float64)
	xout=copy(x);
	taper!(xout,perc)
	return xout
end

function taper!(x, perc=0.0; bperc=perc,eperc=perc)

	nt=size(x,1)
	nttb=min(round(Int,nt*bperc*0.01), nt)
	ntte=min(round(Int,nt*eperc*0.01), nt)
	kb=inv(2.0*round(Int,nt*bperc*0.01)-1)*pi
	ke=inv(2.0*round(Int,nt*eperc*0.01)-1)*pi
	for i in CartesianIndices(size(x))
		if(1≤i[1]≤nttb)
			x[i] *= sin((i[1]-1)*kb)
		end
		if(nt-ntte+1≤i[1]≤nt)
			x[i] *= sin((-i[1]+nt)*ke)
		end
	end
end



"""
Tapering is necessary to be able to input random signal into finite-difference code
Filtering tapering are applied only if the length of the time series is greater than 10
"""
function get_tapered_random_tmax_signal(tgrid; 
					fmin=nothing,
					fmax=nothing,
					tmaxfrac::Float64=1.0,
					dist=Uniform(-2.0, 2.0),
					sparsep=1.0,
					taperperc=20.
					)
	filt_flag=(length(tgrid) > 5) && (!(fmin === nothing)) && (!(fmax===nothing))

	fs = inv(step(tgrid));
	if(filt_flag)
		designmethod = Butterworth(6);
		filtsource = Bandpass(fmin, fmax; fs=fs);
	end

	itind = argmin(abs.(tgrid.-abs(tmaxfrac)*tgrid[end]))
	if(tmaxfrac>0.0)
		its=1:itind
	elseif(tmaxfrac<0.0)
		its=itind:length(tgrid)
	end
	# 20% taper window
	twin = taper(ones(length(its)),taperperc) 
	X = zeros(length(its))
	wavsrc = zeros(length(tgrid)) 
	if(filt_flag) 
		X[:] = rand(dist, length(its)) .* twin
	else
		X[:] = rand(dist, length(its))
	end
	if(sparsep ≠ 1.0)
		Xs=sprandn(length(X), sparsep)
		X[findall(Xs.==0.0)].=0.0
	end
	# band limit
	(filt_flag) && (filt!(X, digitalfilter(filtsource, designmethod), X))
	
	(length(X) ≠ 1) && normalize!(X)
	wavsrc[its] = X
	return wavsrc
end


