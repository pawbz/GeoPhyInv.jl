__precompile__()

module Interferometry

import JuMIT.Grid
import JuMIT.Acquisition
import JuMIT.Data


"""
enhance diffractions in the `TD`

# Keyword Arguments

`λdom::Float64=0.0` : distance between receivers must be greater than twice central wavelength, 2*λdom (Shapiro 2005)
`tlag::Float64=data.tgrid.x[end]-data.tgrid.x[1]` : maximum lag time in the output traces 

"""
function TD_virtual_diff(
			 data::Data.TD;
			 λdom::Float64=0.0,
			 tlag::Float64=data.tgrid.x[end]-data.tgrid.x[1]
			)

	nr = maximum(data.acqgeom.nr);	nss = data.acqgeom.nss;	nt = data.tgrid.nx;
	nfield = data.nfield;
	# normalize the records in data
	datan = Data.TD_normalize(data,:recrms)


	# get unique receiver postions; a virtual source at each
	urpos = Acquisition.Geom_get([data.acqgeom],:urpos)
	nur = Acquisition.Geom_get([data.acqgeom],:nur)

	# central wavelength (use zero for testing)
	println(string("dominant wavelength in the data:\t",λdom))

	rx = Array(Vector{Float64}, nur); rz = Array(Vector{Float64}, nur);
	datmat = zeros(2*nt-1, nur, nur, data.nfield);
	for ifield =1:data.nfield
		# loop over virtual sources
		for irs = 1:nur

			irvec = [];
			# loop over second receiver
			for ir = 1:nur
				rpos = [urpos[1][irs], urpos[2][irs]];
				rpos0 = [urpos[1][ir], urpos[2][ir]];
				δrpos = sqrt((rpos[1] - rpos0[1])^2 + (rpos[2] - rpos0[2])^2)
			
				# the distance between receivers must be greater than 2λdom
				# here λdom is the central wavelength (Shapiro 2005)
				if(δrpos > 2.*λdom)

					# find sources that shoot at these two receivers
					sson = Acquisition.Geom_find(data.acqgeom; rpos=rpos, rpos0=rpos0)
					nsson = count(x->x!=[0],sson);
					if(nsson!=0)
						push!(irvec, ir)
					end
					# stacking over these sources
					for isson=1:length(sson)
						if(sson[isson] != [0])
							datmat[:, ir, irs, ifield] += 
							xcorr(
							  datan.d[isson,ifield][:, sson[isson][2]], 
							  datan.d[isson,ifield][:, sson[isson][1]])
						end
					end
					# normalize depending on the stack
					nsson != 0 ? datmat[:, ir, irs, ifield] /= nsson : nothing
				end
			end
			if(irvec != [])
				rx[irs] = [urpos[2][i] for i in irvec]
				rz[irs] = [urpos[1][i] for i in irvec]
			end
		end
	end
	# virtual source positions 
	sx = [[urpos[2][irs]] for irs in 1:nur];
	sz = [[urpos[1][irs]] for irs in 1:nur];

	# bool array acoording to undef
	ars = [isdefined(rx,irs) ? true : false for irs=1:nur];

	# select only positions that are active
	rx = rx[ars];	rz = rz[ars];
	sx = sx[ars];	sz = sz[ars];
	
	# number of virtual sources
	nvs = length(rx) == length(sx) ? length(sx) : error("some error")

	# geom
	vacqgeom = Acquisition.Geom(sx, sz, rx, rz, nvs, fill(1,nvs), 
			     [length(rx[ir]) for ir=1:length(rx)])

	# tgrid after correlation
	dt = data.tgrid.x[end]-data.tgrid.x[1]
	tgridxcorr=Grid.M1D(-dt, +dt, data.tgrid.δx)
	
	tlag > 0.0 ? tgridcut=Grid.M1D(-tlag, +tlag, data.tgrid.δx) : error("tlag < 0")

	return Data.TD_resamp(
		       Data.TD_urpos(datmat,data.nfield,tgridxcorr,vacqgeom,nur,urpos), 
		       tgridcut)

end



end # module

