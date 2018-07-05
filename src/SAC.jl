
module SAC


using SAC
using Interpolation
using Grid
using ProgressMeter



function process(T)

	npts=maximum(getfield.(T,:npts))

	Tpro=SAC.SACtr[]

	@showprogress 1 "Processing..." for ir in eachindex(T)
		TT=deepcopy(T[ir])
		if(TT.npts == npts)
			# normalize
			SAC.multiply!(TT, inv(vecnorm(TT.t)))

			#SAC.update_headers!(TT)
			push!(Tpro, TT)
		end
	end
	# sort according to distance
	dist=getfield.(Tpro, :gcarc)
	Tpro=Tpro[sortperm(dist)]
	return Tpro
end


function stackamp(T)


	npts=Int64(maximum(getfield.(T,:npts)))
	tgrid=Grid.M1D(Float64(T[1].b), Float64(T[1].e), npts)
	fgrid=Grid.M1D_rfft(tgrid)

	ampstack=zeros(fgrid.nx)

	@showprogress 1 "Processing..." for ir in eachindex(T)
		TT=T[ir]
		if(TT.npts == npts)
			amp=abs.(rfft(normalize(Float64.(TT.t))))

			for i in eachindex(amp)
				ampstack[i]+=amp[i]
			end
		end
	end

	normalize!(ampstack)
	powwav = (ampstack.^2)
	powwavdb = 10. * log10.(powwav./maximum(powwav)) # power in decibel after normalizing
	return powwavdb, fgrid
end

end 
