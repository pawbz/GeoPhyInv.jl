

"""
* `H : number of grid points inside wavelength
	* choose 5 for 4th order scheme
* `epsilon` : Courant number
"""
function check_stability(pa::PFdtd, verbose=true; H=div(20,FD_ORDER), epsilon=inv(sqrt(ndims(pa.c.medium))))
	mgrid=pa.c.medium.mgrid
	freqmin=pa.c.fc[:freqmin]
	freqmax=pa.c.fc[:freqmax]
	attrib_mod=pa.c.attrib_mod
	bounds=pa.c.medium.bounds

	δ=step.(mgrid)
	δt=pa.c.fc[:dt]

	if(isa(attrib_mod, FdtdElastic))
		vmin=minimum(filter(x->x ≠ 0, pa.c.medium[:vs])) # vs condition overrides (other than zero?)
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
