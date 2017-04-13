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


# select virtual source
for irs = 1:


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

