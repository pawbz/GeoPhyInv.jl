module Interferometry

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

	for ifield =1:data.nfield
		dattemp = zeros(2*nt-1, length(urpos), length(urpos));
		sx = [[]];
		# loop over virtual sources
		for irs = 1:length(urpos)
			rx = []
			rz = []
			# loop over second receiver
			for ir = 1:length(urpos)
				# find sources that shoot at these two receivers
				sson = Acquisition.Geom_find((acqgeom; 
					 rpos=[urpos[1][irs], urpos[2][irs]],
					rpos0=[urpos[1][ir], urpos[2][ir]]))
				nsson = length(find(sson));
				if(nss!=0)
					rx[irs] = vcat(rx[irs],urpos[2][ir])
					rz[irs] = vcat(rz[irs],urpos[1][ir])
				end
				# stacking over these sources
				for isson=1:length(sson)
					if(sson[isson])
						dattemp[:, ir, irs] += xcorr(datan[:, ig, isson, ifield],
							datan[:, igs, isson, ifield])
					end
				end
				# normalize depending on the stack
				nsson != 0 ? dattemp[:, ir, irs] /= nsson : nothing
			end
		end
	end


	dattemp = zeros(nt, ng, ng, ifield) # each receiver is turned into a virtual source
	for ifield = 1:data.nfield, is = 1:nss, ig = 1:nr
		dattemp[:, ig, is, ifield] += xcorr(datan[:, ig, is, ifield],
							datan[:, ig, is, ifield])
	end
	acqvirtual = data.acqgeom;
	acqvirtual.nss = maximum(acqvirtual.nr);
	for 
	acqvirtual.sx = reshape(sx, 1, acqvirtual.nss);
	szall = reshape(sz, 1, nss);

	dataout = TD(dattemp,data.nfield,data.tgrid,acqvirtual)
	return Data.TD()

end



end # module

