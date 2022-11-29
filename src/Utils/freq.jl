
#=

"""
Return the normalized and stacked (along dim=1) power spectra of time series x in dB 
as an AxisArray.
"""
function stacked_power_spectra(x, fs=1.0)
	# allocate rfft domain
	cx=rfft(x,[1]); # rfft along first dimension
	
	abs2.(cx) #  
	for i in CartesianIndices(size(cx))
		rmul!()
		if(1≤i[1]≤nttb)
			x[i] *= sin((i[1]-1)*kb)
		end
		if(nt-ntte+1≤i[1]≤nt)
			x[i] *= sin((-i[1]+nt)*ke)
		end
	end


end
=#


function findfreq(
		  x::AbstractArray{T, ND} where {T<:AbstractFloat},
		  tgrid;
		  attrib::Symbol=:peak,
		  threshold::Float64=-50.
		  ) where {ND}

	cx=rfft(x,[1]);
	fgrid=FFTW.rfftfreq(length(tgrid),inv(step(tgrid))) # corresponding fgrid

	ax = (abs.(cx).^2); # power spectrum in dB

	if(maximum(ax) == 0.0)
		@warn("x is zero"); return 0.0
	else 
		ax /= maximum(ax);
		ax = 10. .* log10.(ax)
	end

	if(attrib == :max)
		iii=findlast(ax .>= threshold)
	elseif(attrib == :min)
		iii=findfirst(ax .>= threshold)
	elseif(attrib == :peak)
		iii=argmax(ax)
	end
	ii=CartesianIndices(size(ax))[iii][1]
	return fgrid[ii]

end



function stacked_spectrum(dobs; fs=1.0)
	D=rfft(dobs,1) # perform rfft
	window=zeros(size(D,1))        
	nr=size(dobs,2)        
	# stack the spectrum of data        
	for j in 1:nr                
		for i in eachindex(window)
			window[i] += (abs2(D[i,j]))
		end
	end        
	rmul!(window,inv(maximum(abs,window)))
	window = 10. * log10.(window)        
	return Float64.(FFTW.rfftfreq((size(dobs,1)),fs)), window
end


