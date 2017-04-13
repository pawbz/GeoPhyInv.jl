module Interferometry

import SIT.Acquisition
import SIT.Data


"""
enhance diffractions in the `TD`
"""
function TD_virtual_diff(
			 data::Data.TD,
			)

	nr = maximum(data.acqgeom.nr);	nss = data.acqgeom.nss;	nt = data.tgrid.nx;
	nfield = data.nfield;
	# normalize the records in data
	datan = Data.TD_normalize(data,:recrms)


	# get unique receiver postions; a virtual source at each
	urpos = Acquisition.Geom_get(data.acqgeom,:urpos)
	nur = Acquisition.Geom_get(data.acqgeom,:nur)

	rx = Array(Vector{Float64}, nur); rz = Array(Vector{Float64}, nur);
	for ifield =1:data.nfield
		dattemp = zeros(2*nt-1, nur, nur);
		# loop over virtual sources
		for irs = 1:nur

			irvec = [];
			# loop over second receiver
			for ir = 1:nur
				# find sources that shoot at these two receivers
				sson = Acquisition.Geom_find(data.acqgeom; 
					 rpos=[urpos[1][irs], urpos[2][irs]],
					rpos0=[urpos[1][ir], urpos[2][ir]])
				nsson = count(x->x!=[0],sson);
				if(nss!=0)
					push!(irvec, ir)
				end
				# stacking over these sources
				for isson=1:length(sson)
					if(sson[isson] != [0])
						dattemp[:, ir, irs] += 
							xcorr(datan.d[:, sson[isson][2], isson, ifield],
						   datan.d[:, sson[isson][1], isson, ifield])
					end
				end
				# normalize depending on the stack
				nsson != 0 ? dattemp[:, ir, irs] /= nsson : nothing

			end
			if(irvec != [])
				println(irvec)
				rx[irs] = [urpos[2][i] for i in irvec]
				rz[irs] = [urpos[1][i] for i in irvec]
			end

		end
	end

	println(rx, typeof(rx))
	sx = Array(Vector{Float64}, nur); sz = Array(Vector{Float64}, nur);
	return rx, rz

end



end # module

