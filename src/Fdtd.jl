module Fdtd

import SIT.Grid
import SIT.Models
import SIT.F90libs
import SIT.Acquisition
import SIT.Data
import SIT.Gallery
import SIT.DSP


"""
test function
"""
function __init__()

end

"""
born modeling

# Arguments
* `tgridmod`: time-domain grid for modeling
* `tgrid`: time grid for output data 
"""
function fdtd_born_mod(;
		  jobname::AbstractString = "Hello",
		  npropwav::Int64 = 1, 
		  model::Models.Seismic = Gallery.Seismic(:seismic_homo1),
		  tgridmod::Grid.M1D = Gallery.M1D(:seismic_homo1),
		  tgrid::Grid.M1D = tgridmod,
		  acqgeom::Acquisition.Geom = Gallery.Geom(:seismic_homo1),
		  acqsrc::Acquisition.Src = Gallery.Src(:seismic_homo1),
		  src_nsmul::Int64 = 1,
		  src_flags::AbstractString = "[BILINEAR]",
		  recv_flags = "[BILINEAR]"
		 )

#if(tgrid.nx .lt. src_nt) call abort_msg("fd2_mod: tgrid.nx .lt. src_nt")

#! no modeling if source wavelet is zero
#if(maxval(abs(src_wavelets)) .lt. tiny(rzero_de)) then
#        return
#endif
tgridmod.δx > tgrid.δx ? error("output time grid sampling finer than modeling") :


# find maximum and minimum frequencies in the source wavelets
freqmin = DSP.findfreq(acqsrc.wav[:,1,1],acqsrc.tgrid,attrib=:min) 
freqmax = DSP.findfreq(acqsrc.wav[:,1,1],acqsrc.tgrid,attrib=:max) 

# check spatial sampling
ds_temp=Models.χ([minimum(model.χvp)],model.vp0,-1)[1]/5.0/freqmax;
ds_max = maximum([model.mgrid.δx, model.mgrid.δz])
all(ds_max .> ds_temp) ? 
		error(string("spatial sampling\t",ds_max,"\ndecrease spatial sampling below:\t",ds_temp)) :
		println("spatial sampling can be as high as:\t",ds_temp)

# check time sampling
ds_min = minimum([model.mgrid.δx, model.mgrid.δz])
dt_temp=0.5*ds_min/Models.χ([(maximum(model.χvp))],model.vp0,-1)[1]
all(tgridmod.δx .> dt_temp) ? 
		error(string("time sampling\t",tgridmod.δx,"\ndecrease time sampling below:\t",dt_temp)) :


recv_n = maximum(acqgeom.nr);
src_nseq = acqgeom.nss

# some allocations
recv_out = zeros(tgridmod.nx*recv_n*npropwav*src_nseq)
snaps_in = zeros(tgridmod.nx)
snaps_out = ones(tgridmod.nx)
grad_modtt = zeros(model.mgrid.nz, model.mgrid.nx)
grad_modrr = zeros(model.mgrid.nz, model.mgrid.nx)


# save Greens function in the background velocity model 
#acqgeomb = Acquisition.Geom(acqgeom, mgrid,:recborder) 
npropwav = 1;
records = fdtd_mod(jobname="for",model=model,tgridmod=tgridmod,
		   acqgeom=[acqgeom], acqsrc=[acqsrc],
	npropwav=npropwav, 
	)

# receivers to sources

#datamod = Data.TD(reshape(records,tgridmod.nx,recv_n*npropwav,src_nseq),tgridmod,acqgeom)

return reshape(records,tgridmod.nx,recv_n*npropwav,src_nseq)

#return Data.TD_resamp(datamod, tgrid)

#isapprox(maximum(abs(recv_out)),0.0) && warn("receiver records are zeros")


#return records
end



"""
As forward modeling method, the 
finite-difference method is employed. 
It uses a discrete version of the two-dimensional isotropic acoustic wave equation.
As shown in  
Figure~ef{fdshmesh_acou}, a staggered-grid mesh is used 
# Arguments
* `jobname` : dominant frequency
* `npropwav` : number of wavefields propagating independently in the same medium time-domain grid
* `model` : Seismic Model
* `model0` : Background Seismic Model used for Born modeling 
* `tgridmod` : time grid for modeling
* `tgrid` : time grid for data output
* `prop_flags`: flags that combine propagating wavefields
"""
function fdtd_mod(;
		  jobname::AbstractString = "Hello",
		  npropwav::Int64 = 1, 
		  model::Models.Seismic = Gallery.Seismic(:seismic_homo1),
		  model0::Models.Seismic = Gallery.Seismic(:seismic_homo1),
		  tgridmod::Grid.M1D = Gallery.M1D(:seismic_homo1),
		  tgrid::Grid.M1D = tgridmod,
		  acqgeom::Array{Acquisition.Geom} = [Gallery.Geom(:seismic_homo1)],
		  acqsrc::Array{Acquisition.Src} = [Gallery.Src(:seismic_homo1)],
		  src_flags::AbstractString = "[BILINEAR]",
		  recv_flags = "[BILINEAR]",
		  prop_flags = "[BILINEAR]",
		  verbose::Bool = false
		 )

maximum(tgridmod.x) < maximum(acqsrc[1].tgrid.x) ? error("modeling time is less than source time") :

tgridmod.δx > tgrid.δx ? error("output time grid sampling finer than modeling") :

# find maximum and minimum frequencies in the source wavelets
freqmin = DSP.findfreq(acqsrc[1].wav[:,1,1],acqsrc[1].tgrid,attrib=:min) 
freqmax = DSP.findfreq(acqsrc[1].wav[:,1,1],acqsrc[1].tgrid,attrib=:max) 

# check spatial sampling
ds_temp=Models.χ([minimum(model.χvp)],model.vp0,-1)[1]/5.0/freqmax;
ds_max = maximum([model.mgrid.δx, model.mgrid.δz])
all(ds_max .> ds_temp) ? 
		error(string("spatial sampling\t",ds_max,"\ndecrease spatial sampling below:\t",ds_temp)) :
		verbose ? println("spatial sampling can be as high as:\t",ds_temp) : nothing 

# check time sampling
ds_min = minimum([model.mgrid.δx, model.mgrid.δz])
dt_temp=0.5*ds_min/Models.χ([(maximum(model.χvp))],model.vp0,-1)[1]
all(tgridmod.δx .> dt_temp) ? 
		error(string("time sampling\t",tgridmod.δx,"\ndecrease time sampling below:\t",dt_temp)) :
		verbose ? println("time sampling can be as high as:\t",dt_temp) : nothing



#! no modeling if source wavelet is zero
#if(maxval(abs(src_wavelets)) .lt. tiny(rzero_de)) then
#        return
#endif

# all the propagating wavefield should have same sources, receivers, ? check that?

# check dimension of model and model0

# check if all sources are receivers are inside model

length(acqgeom) != npropwav ? error("acqgeom size") : nothing
length(acqsrc) != npropwav ? error("acqsrc size") : nothing
any([getfield(acqgeom[ip],:nss) != getfield(acqsrc[ip],:nss) for ip=1:npropwav])  ? error("different supersources") : nothing
any([getfield(acqgeom[ip],:ns) != getfield(acqsrc[ip],:ns) for ip=1:npropwav])  ? error("different sources") : nothing

recv_n = maximum(acqgeom[1].nr);
src_nsmul = maximum(acqgeom[1].ns);
src_nseq = acqgeom[1].nss;

if(verbose)
	println(string("number of receivers:\t",recv_n))	
	println(string("number of sources:\t",src_nsmul))	
	println(string("number of super sources:\t",src_nseq))	
end

recv_out = zeros(tgridmod.nx*recv_n*npropwav*src_nseq)
snaps_in = zeros(tgridmod.nx) # dummy
snaps_out = ones(tgridmod.nx) # dummy
grad_modtt = zeros(model.mgrid.nz, model.mgrid.nx) # dummy
grad_modrr = zeros(model.mgrid.nz, model.mgrid.nx) # dummy

# extend models in the PML layers
exmodel = Models.Seismic_extend(model);
exmodel0 = Models.Seismic_extend(model0);


ccall( (:fdtd_mod, F90libs.fdtd), Void,
      (Ptr{UInt8}, Ref{Int64},       Ptr{Float64}, Ptr{Float64},
       Ptr{Float64}, Ptr{Float64},
       Ref{Int64}, Ref{Int64},       Ref{Float64}, Ref{Float64},
       Ref{Float64}, Ref{Float64},   Ref{Int64}, Ptr{UInt8},
       Ref{Int64}, Ref{Float64},     Ref{Int64}, Ptr{Float64},
       Ref{Int64}, Ref{Int64},       Ptr{Float64}, Ptr{Float64},
       Ptr{UInt8}, Ref{Int64},       Ptr{Float64}, Ptr{Float64},
       Ptr{Float64}, Ptr{UInt8}, Ptr{UInt8},
       Ptr{Float64}, Ptr{Float64},
       Ptr{Float64}, Ptr{Float64}
       ),
      jobname, npropwav,     Models.χ(exmodel.χKI, exmodel.K0I,-1), Models.χ(exmodel.χρI, exmodel.ρ0I,-1),
      Models.χ(exmodel0.χKI, exmodel0.K0I,-1), Models.χ(exmodel0.χρI, exmodel0.ρ0I,-1),
      model.mgrid.nx, model.mgrid.nz,    model.mgrid.δx, model.mgrid.δz,
      model.mgrid.x, model.mgrid.z,      model.mgrid.npml-1, # reduce npml by one, see fdtd.f90
      "model.mgrid_abs_trbl",
      tgridmod.nx, tgridmod.δx,    acqsrc[1].tgrid.nx, Acquisition.get_vecSrc(:wav,acqsrc),
      src_nseq, src_nsmul, Acquisition.get_vecGeom(:sx,acqgeom), Acquisition.get_vecGeom(:sz,acqgeom),
      src_flags, recv_n, Acquisition.get_vecGeom(:rx,acqgeom), Acquisition.get_vecGeom(:rz,acqgeom),
      recv_out, recv_flags,  prop_flags,
      snaps_in, snaps_out,
      grad_modtt, grad_modrr,
     )

# check if ccall return zeros
isapprox(maximum(abs(recv_out)),0.0) && warn("recv_out are zeros")

# return after forming a vector and resampling
nd = src_nseq*tgridmod.nx*recv_n; # number of samples for each iprop
return [Data.TD_resamp(Data.TD(reshape(recv_out[1+(iprop-1)*nd : iprop*nd],tgridmod.nx,recv_n,src_nseq),
		       tgridmod, acqgeom[1]), tgrid) for iprop in 1:npropwav]

# return without resampling for testing
#return [Data.TD(reshape(recv_out[1+(iprop-1)*nd : iprop*nd],tgridmod.nx,recv_n,src_nseq),
#		       tgridmod, acqgeom[1]) for iprop in 1:npropwav]
end

end # module
