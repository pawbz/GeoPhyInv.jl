__precompile__()

"""
This module defines the data types related to seismic data:
* `TD` : time domain representation

It also provides methods that apply source and receiver filters onto 
seismic data.
"""
module Data

import JuMIT.Acquisition
import JuMIT.Grid
import JuMIT.Interpolation
import JuMIT.Coupling
import JuMIT.DSP

"""
Time domain representation of Seismic Data.

# Fields

* `d::Array{Array{Float64,2},2}` : data 
* `fields::Vector{Symbol}` : components recorded at each receiver
* `tgrid::Grid.M1D` : grid to represent time
* `acqgeom::Acquisition.Geom` : acquisition geometry used to generate the data
"""
type TD
	d::Array{Array{Float64,2},2}
	fields::Vector{Symbol}
	tgrid::Grid.M1D
	acqgeom::Acquisition.Geom
	"adding conditions that are to be false while construction"
	TD(d, fields, tgrid, acqgeom) = 
		any([
		  any([fields[ifield] ∉ [:P, :Vx, :Vz] for ifield in 1:length(fields)]),
		  length(fields) == 0,
		  broadcast(size,d) != [(tgrid.nx,acqgeom.nr[iss]) for iss=1:acqgeom.nss, ifield=1:length(fields)]
		  ]) ? 
		error("error in TD construction") : new(d, fields, tgrid, acqgeom)

end


"""
Method to resample data in time.

# Arguments

* `data::TD` : input data of type `TD`
* `tgrid::Grid.M1D` : resampling in time according to this time grid

# Return

* data after resampling as `TD`
"""
function TD_resamp(data::TD,
		tgrid::Grid.M1D
		)
	nss = data.acqgeom.nss
	nr = data.acqgeom.nr
	dataout = TD(
	      [zeros(tgrid.nx,data.acqgeom.nr[iss]) for iss=1:nss, ifield=1:length(data.fields)],
	      length(data.fields),tgrid,data.acqgeom)
	TD_resamp!(dataout, data)
	return dataout
end


"""
Method to resample data in time.

# Arguments

* `data::TD` : input data of type `TD`
* `dataout::TD` : preallocated data of type `TD` that is modified
"""
function TD_resamp!(dataout::TD, data::TD)
	# check if datasets are similar
	nss = data.acqgeom.nss
	nr = data.acqgeom.nr
	for ifield = 1:length(data.fields), iss = 1:nss, ir = 1:nr[iss]
		din = data.d[iss, ifield][:, ir]
		xin = data.tgrid.x
		xout = dataout.tgrid.x
		dout = similar(xout)
		Interpolation.interp_spray!(xin, din, xout, dout, :interp, :B1 )

		dataout.d[iss, ifield][:,ir] = copy(dout)
	end
	return dataout
end


"""
Method used to preallocate `TD` with zeros.

# Arguments

* `fields::Vector{Symbol}` : number of components
* `tgrid::Grid.M1D` : time domain grid
* `acqgeom::Acquisition.Geom` : acquisition geometry

# Return

* data with zeros as `TD`
"""
function TD_zeros(fields::Vector{Symbol}, tgrid::Grid.M1D, acqgeom::Acquisition.Geom)
	return TD([zeros(tgrid.nx,acqgeom.nr[iss]) for iss=1:acqgeom.nss, ifield=1:length(fields)],fields,
	   deepcopy(tgrid),deepcopy(acqgeom)) 
end
"Same as `TD_zeros`, except for returning ones"
function TD_ones(fields::Vector{Symbol}, tgrid::Grid.M1D, acqgeom::Acquisition.Geom) 
	return TD([ones(tgrid.nx,acqgeom.nr[iss]) for iss=1:acqgeom.nss, ifield=1:length(fields)],
	   fields,deepcopy(tgrid),deepcopy(acqgeom)) 
end


"""
Returns bool depending on if input `data::TD` has all zeros or not.
"""
function TD_iszero(data::TD)
	return maximum(broadcast(maximum,abs,data.d)) == 0.0 ? true : false
end

"""
Time reverse the records of each receiver in `TD` 

# Arguments

* `data::TD` : input data that is modified
"""
function TD_tr!(data::TD)
	data.d = copy([flipdim(data.d[i,j],1) for i in 1:data.acqgeom.nss, j in 1:length(data.fields)]);
end

"""
Returns dot product of data.

# Arguments 

* `data1::TD` : data 1
* `data2::TD` : data 2

# Return

* dot product as `Float64`
"""
function TD_dot(data1::TD, data2::TD)
	dotd = 0.0;
	for ifield = 1:length(data1.fields), iss = 1:data1.acqgeom.nss, ir = 1:data1.acqgeom.nr[iss]
		dotd += dot(data1.d[iss, ifield][:, ir],data2.d[iss, ifield][:, ir])
	end
	return dotd
end


"""
Normalize time-domain seismic data.

# Arguments 

* `data::TD` : input data
* `attrib::Symbol` : decide kind of normalization
  * `=:recrms` the record at every receiver is normalized with its RMS value
  * `=:recmax` the record at every receiver is normalized with its maximum value

# Return

* normalized data as `TD`
"""
function TD_normalize(data::TD, attrib::Symbol)
	nr = data.acqgeom.nr;	nss = data.acqgeom.nss;	nt = data.tgrid.nx;
	datan = deepcopy(data);
	for ifield = 1:length(data.fields), iss = 1:nss, ir = 1:nr[iss]
		if(attrib == :recrms)
			nval = sqrt(mean(datan.d[iss, ifield][:,ir].^2.))
		elseif(attrib == :recmax)
			nval = maximum(datan.d[iss, ifield][:,ir])
		else
			error("invalid attrib")
		end

		# normalize
		datan.d[iss, ifield][:, ir] = 
		isequal(nval, 0.0) ? zeros(nt) : datan.d[iss, ifield][:, ir]./nval  
	end
	return datan
end

"""
Construct TD using data at all the unique receiver positions
for all supersources.

* `d::Array{Float64}` : the data matrix ordered in order such that time-domain modelling schemes are fast, i.e., [irec,ifield,it,nss]

"""
function TD_urpos(d::Array{Float64}, 
		   fields::Vector{Symbol}, 
		   tgrid::Grid.M1D, 
		   acq::Acquisition.Geom,
		   nur::Int64,
		   urpos::Tuple{Array{Float64,1},Array{Float64,1}
		  }
		   )
	dout = [zeros(tgrid.nx,acq.nr[iss]) for iss=1:acq.nss, ifield=1:length(fields)] 

	for ifield=1:length(fields), iss=1:acq.nss, ir=1:acq.nr[iss]
		# find index in urpos
		irr=find([[urpos[1][i]-acq.rz[iss][ir],
		       urpos[2][i]-acq.rx[iss][ir]] == [0., 0.,] for i in 1:nur])

		dout[iss, ifield][:,ir] = d[irr[1],ifield, :,iss] 
	end

	return TD(dout, fields, tgrid, acq)

end

"""
Apply source and receiver coupling functions to TD.
Currently, only source filters are applied.

# Arguments

* `s::TD` : input data
* `r::TD` : input data
* `w::Coupling.TD` : input source and receiver filters
* `attrib::Symbol` : attribute to 
  * `=:s` to apply `w` to `r` and modify `s`
  * `=:r` to apply adjoint of `w` to `s` and modify `r`
  * `=:w` modify `w` using `r` and `s`

TODO: need to work on parallelization and speed up here
"""
function TDcoup!(
               s::TD,
	       r::TD,
	       w::Coupling.TD,
	       attrib::Symbol
	       )
	nr = r.acqgeom.nr;	nss = r.acqgeom.nss;	nt = r.tgrid.nx;
	fields = (w.fields == r.fields == s.fields) ? w.fields : error("different fields")
	sv=zeros(s.tgrid.nx)
	rv=zeros(r.tgrid.nx)
	wv=zeros(w.tgridssf.nx)
	for ifield = 1:length(fields), iss = 1:nss, ir = 1:nr[iss]
		# receiver coupling
	#	DSP.fast_filt!(s.d[iss, ifield][:, ir],r.d[iss, ifield][:, ir],
#		 w.rf[iss, ifield][:,ir],:s)
		# source coupling
		sv=s.d[iss, ifield][:, ir]
		rv=r.d[iss, ifield][:, ir]
		wv=w.ssf[iss, ifield]
		DSP.fast_filt!(sv, rv, wv, attrib)
		s.d[iss, ifield][:, ir]=copy(sv)
		r.d[iss, ifield][:, ir]=copy(rv)
		w.ssf[iss, ifield][:]=copy(wv)
	end
end # TDcoup


"""
Apply different weighting functions to `TD`. Use this method to create a 
data preconditioning matrix

# Arguments Modified

* `dw::TD` : 

# Keyword Arguments

* `offsetlim::Vector{Float64}=[-Inf,Inf]` : [xoffsetlim, zoffsetlim], where the records with offsets > `offsetlim` are given zero weight 
* `tlim::Vector{Float64}=[dw.tgrid.x[1], dw.tgrid.x[end]]` : [tminimum, tmaximum], time mute window 
* `offsetpow::Vector{Float64}=[0.0,0.0]` : 
* `tpow::Float64=0.0` :
* `ttaperperc::Float64=0.` : taper window percentage for time window

* NOTE: if more than one simultaneous source are present, their mean position is considered to calculate offset.
"""
function TD_weight!(
		    dw::TD;
		    offsetlim::Vector{Float64}=[-Inf,Inf],
		    tlim::Vector{Float64}=[dw.tgrid.x[1], dw.tgrid.x[end]],
		    offsetpow::Vector{Float64}=[0.0,0.0],
		    tpow::Float64=0.0,
		    ttaperperc::Float64=0.,
		  )

	tvecexp = dw.tgrid.x - dw.tgrid.x[1]
	tmaxI = maximum(abs(tvecexp))^(-1)
	nt = dw.tgrid.nx
	fields=dw.fields
	nss=dw.acqgeom.nss
	rx=dw.acqgeom.rx; rz=dw.acqgeom.rz
	sx=dw.acqgeom.sx; sz=dw.acqgeom.sz
	nr=dw.acqgeom.nr; ns=dw.acqgeom.ns

	itlim = sort(broadcast(indmin,[abs(dw.tgrid.x-tlim[i]) for i in 1:2]))
	twin=zeros(nt)
	twin[itlim[1] : itlim[2]] = DSP.taper(ones(itlim[2]-itlim[1]+1),ttaperperc) 

	for ifield = 1:length(fields), iss = 1:nss
		zo = sqrt((rz[iss][:]-mean(sz[iss])).^2) # offsets computed using mean of source position
		xo = sqrt((rx[iss][:]-mean(sx[iss])).^2)
		zomaxI = maximum(zo)^(-1)
		xomaxI = maximum(xo)^(-1)
		for ir = 1:nr[iss]
			inoffsetlim = ( (abs(zo[ir]) < abs(offsetlim[2])) & (abs(xo[ir]) < abs(offsetlim[1])) )
			if(inoffsetlim)
				# apply xoffset weighting 
				if(!isinf(xomaxI))
					dw.d[iss,ifield][:,ir] .*= exp(xo[ir]*offsetpow[1]*xomaxI) 
				end
				# apply zoffset weighting
				if(!isinf(zomaxI))
					dw.d[iss,ifield][:,ir] .*= exp(zo[ir]*offsetpow[2]*zomaxI)
				end
			else
				dw.d[iss,ifield][:,ir] = 0.
			end
			for it=1:nt
				dw.d[iss,ifield][it,ir] *= twin[it] * exp(tvecexp[it]*tpow*tmaxI)
			end
		end
	end

end

end # module
