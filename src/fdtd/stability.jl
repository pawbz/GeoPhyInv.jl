

"""
* `H : number of grid points inside wavelength
	* choose 5 for 4th order scheme?
* `epsilon` : Courant number
"""
function check_fd_stability(bounds, mgrid, tgrid, freqmin, freqmax, attrib_mod, verbose=true, H=5, epsilon=0.5)

	δ=step.(mgrid)
	δt=step(tgrid)

	if(isa(attrib_mod, FdtdElastic))
		vmin=bounds[:vs][1] # vs condition overrides 
		vmax=sqrt(bounds[:vp][2]^2 + bounds[:vs][2]^2) # see Virieux (1986)

	else
		vmin=bounds[:vp][1]
		vmax=bounds[:vp][2]
	end

	# check spatial sampling, i.e., number of grid points in minimum wavelength
	if(isa(attrib_mod, FdtdElastic))
		δs_temp=round(vmin/H/freqmax,digits=2); 
	else
		δs_temp=round(vmin/H/freqmax,digits=2);
	end
	δs_max = maximum(δ)
	str1=@sprintf("%0.2e",δs_max)
	str2=@sprintf("%0.2e",δs_temp)
	if(str1 ≠ str2)
		if(all(δs_max .> δs_temp)) 
			@warn "decrease maximum spatial sampling ($str1) below $str2"
		else
			if(verbose)
				@info "spatial sampling ($str1) can be as high as $str2"
			end
		end
	end

	# check time sampling, i.e., Courant criteria 
	δs_min = minimum(δ)
	δt_temp=epsilon*δs_min/vmax
	str1=@sprintf("%0.2e", δt)
	str2=@sprintf("%0.2e", δt_temp)
	if(str1 ≠ str2)
		if(all(δt .> δt_temp))
			@warn "decrease time sampling ($str1) below $str2"
		else
			if(verbose)
				@info "time sampling ($str1) can be as high as $str2"
			end
		end
	end
end
