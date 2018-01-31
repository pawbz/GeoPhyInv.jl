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
import JuMIT.Conv
import JuMIT.Interpolation
import JuMIT.Coupling
import JuMIT.Misfits
import JuMIT.DSP
using DSP

"""
Time domain representation of Seismic Data.

# Fields

* `d::Array{Array{Float64,2},2}` : data 
* `fields::Vector{Symbol}` : components recorded at each receiver
* `tgrid::Grid.M1D` : grid to represent time
* `acqgeom::Acquisition.Geom` : acquisition geometry used to generate the data
"""
type TD
	d::Matrix{Matrix{Float64}}
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

"Compare if two `TD`'s  are equal"
function Base.isequal(dat1::TD, dat2::TD)
	fnames = fieldnames(TD)
	vec=[(isequal(getfield(dat1, name),getfield(dat2, name))) for name in fnames]
	return all(vec)
end

"""
Return if two `TD`'s have same dimensions and bounds.
"""
function Base.isapprox(dat1::TD, dat2::TD)
	vec=([ 
		isequal(dat1.tgrid, dat2.tgrid),
		isequal(dat1.fields, dat2.fields),
		isequal(dat1.acqgeom, dat2.acqgeom, :receivers), # only receivers have to be the same
       		(size(dat1.d)==size(dat2.d)), 
		])
	vec2=[size(dat1.d[iss,ifield])==size(dat2.d[iss,ifield]) for iss=1:dat1.acqgeom.nss, ifield=1:length(dat1.fields)]
	return (all(vec) & all(vec2))
end



"""
Copy `TD`'s, which are similar.
"""
function Base.copy!(dataout::TD, data::TD)
	if(isapprox(dataout, data))
		dout=getfield(dataout, :d)
		din=getfield(data, :d)
		for iss=1:data.acqgeom.nss, ifield=1:length(data.fields)
			ddout=dout[iss,ifield]
			ddin=din[iss,ifield]
			copy!(ddout,ddin)
		end
		return dataout
	else
		error("attempt to copy dissimilar data")
	end
end



"""
Method to resample data in time.

# Arguments

* `data::TD` : input data of type `TD`
* `tgrid::Grid.M1D` : resampling in time according to this time grid

# Return

* data after resampling as `TD`
"""
function interp(data::TD,
		tgrid::Grid.M1D
		)
	nss = data.acqgeom.nss
	nr = data.acqgeom.nr
	dataout = TD(
	      [zeros(tgrid.nx,data.acqgeom.nr[iss]) for iss=1:nss, ifield=1:length(data.fields)],
	      data.fields,tgrid,data.acqgeom)
	interp_spray!(dataout, data)
	return dataout
end


"""
Method to resample data in time.
Can reduce allocations =========

# Arguments

* `data::TD` : input data of type `TD`
* `dataout::TD` : preallocated data of type `TD` that is modified
"""
function interp_spray!(dataout::TD, data::TD, attrib=:interp, Battrib=:B1)
	# check if datasets are similar
	nss = data.acqgeom.nss
	nr = data.acqgeom.nr
	xin=data.tgrid.x
	xout=dataout.tgrid.x
	for ifield = 1:length(data.fields), iss = 1:nss
		dat=data.d[iss,ifield]
		dato=dataout.d[iss,ifield]
		for ir = 1:nr[iss]
			din=view(dat,:,ir)
			dout=view(dato,:,ir)
	 		Interpolation.interp_spray!(xin, din, xout, dout,attrib,Battrib)
		end
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
function Base.fill!(data::TD, k::Float64)
	for iss=1:data.acqgeom.nss, ifield=1:length(data.fields)
		data.d[iss,ifield][:]=k 
	end
end
function TD_zeros(d::TD)
	return TD([zeros(d.tgrid.nx,d.acqgeom.nr[iss]) for iss=1:d.acqgeom.nss, ifield=1:length(d.fields)],d.fields,
	   deepcopy(d.tgrid),deepcopy(d.acqgeom)) 
end


"Same as `TD_zeros`, except for returning ones"
function TD_ones(fields::Vector{Symbol}, tgrid::Grid.M1D, acqgeom::Acquisition.Geom) 
	return TD([ones(tgrid.nx,acqgeom.nr[iss]) for iss=1:acqgeom.nss, ifield=1:length(fields)],
	   fields,deepcopy(tgrid),deepcopy(acqgeom)) 
end


"""
Returns bool depending on if input `data::TD` has all zeros or not.
"""
function Base.iszero(data::TD)
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
function Base.dot(data1::TD, data2::TD)
	if(isapprox(data1, data2))
		dotd = 0.0;
		for ifield = 1:length(data1.fields), iss = 1:data1.acqgeom.nss 
			for ir = 1:data1.acqgeom.nr[iss], it = 1:data1.tgrid.nx
				dotd += data1.d[iss, ifield][it, ir] * data2.d[iss, ifield][it, ir]
			end
		end
		return dotd
	else
		error("cannot dot dissimilar datasets")
	end
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
function TD_normalize(data::TD, attrib::Symbol=:recrms)
	datan = deepcopy(data);
	TD_normalize!(datan, attrib)
	return datan
end
function TD_normalize!(data::TD, attrib::Symbol=:recrms)
	nr = data.acqgeom.nr;	nss = data.acqgeom.nss;	nt = data.tgrid.nx;
	for ifield = 1:length(data.fields), iss = 1:nss
		dd=data.d[iss, ifield]
		scs=vecnorm(dd,2)

		for ir = 1:nr[iss]
			ddv=view(dd, :, ir)

			if(attrib == :recrms)
				sc=vecnorm(ddv,2)
				scale!(ddv,inv(sc))
			elseif(attrib == :recmax)
				sc=vecnorm(ddv,Inf)
				scale!(ddv,inv(sc))
			elseif(attrib == :srcrms)
				scale!(ddv,inv(scs))
			else
				error("invalid attrib")
			end
		end
	end
	return data
end


function TD_filter!(data::TD; fmin=nothing, fmax=nothing)

	if((fmin===nothing) && (fmax===nothing))
		return nothing
	end
	fs = 1/ data.tgrid.δx;
	designmethod = Butterworth(4);
	filtsource = Bandpass(fmin, fmax; fs=fs);

	nr = data.acqgeom.nr;	nss = data.acqgeom.nss;	nt = data.tgrid.nx;
	for ifield = 1:length(data.fields), iss = 1:nss
		dd=data.d[iss, ifield]
		scs=vecnorm(dd,2)

		for ir = 1:nr[iss]
			ddv=view(dd, :, ir)

			filt!(ddv, digitalfilter(filtsource, designmethod), ddv)
		end
	end
	return data
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



"""
Calculate the distance between the observed data `y` and the calculated data `x`.
The time grid of the observed data can be different from that of the modelled data.
The acqistion geometry of both the data sets should be the same.

If `J` is the distance, the gradient of the misfit w.r.t to the calculated data is returned as `dJx`
* `w` used for data preconditioning
* `coupling` source and receiver coupling functions
"""
mutable struct Param_error
	x::TD # modelled data
	y::TD # observed data
	w::TD # weights
	xr::TD # modelled data after resampling
	dJxr::TD # gradient w.r.t. xr
	dJx::TD # gradient w.r.t. x
	dJssf::Matrix{Vector{Float64}} # gradient w.r.t source filters
	ynorm::Float64 # normalize functional with energy of y
	coup::Coupling.TD
	paconvssf::Conv.Param{Float64,1} # convolutional model for source filters
	paconvrf::Conv.Param{Float64,1} # convolutional model for receiver filters
	func::Function # function to compute misfit
	dJxc::Vector{Float64} # temp array of length nt
end

function Param_error(x,y,w=nothing, coup=nothing)

	if(coup===nothing)
		 coup=Coupling.TD_delta(y.tgrid,[0.1,0.1],[0.0, 0.0], [:P], y.acqgeom)
	end

	paconvssf=Conv.Param(ntwav=coup.tgridssf.nx, 
		      ntd=y.tgrid.nx, 
		      ntgf=y.tgrid.nx, 
		      wavlags=coup.ssflags, 
		      dlags=[y.tgrid.nx-1, 0], 
		      gflags=[y.tgrid.nx-1, 0])

	paconvrf=Conv.Param(ntwav=coup.tgridrf.nx, 
		      ntd=y.tgrid.nx, 
		      ntgf=y.tgrid.nx, 
		      wavlags=coup.rflags, 
		      dlags=[y.tgrid.nx-1, 0], 
		      gflags=[y.tgrid.nx-1, 0])

	dJssf=deepcopy(ssf)

	if(w===nothing) 
		w=Data.TD_ones(y.fields,y.tgrid,y.acqgeom)
	end
	!(isapprox(w,y)) && error("weights have to be similar to y")

	!(isequal(x.acqgeom, y.acqgeom)) && error("observed and modelled data should have same acqgeom")
	!(isequal(x.fields, y.fields)) && error("observed and modelled data should have same fields")

	xr=TD_zeros(y)
	dJxr=TD_zeros(y)
	dJx=TD_zeros(x)

	ynorm = dot(y, y)
	(ynorm == 0.0) && error("y cannot be zero")

	func=Misfits.error_squared_euclidean!
	dJxc=zeros(y.tgrid.nx)

	return Param_error(x,y,w,xr,dJxr,dJx,dJssf,ynorm,coup,paconvssf, paconvrf, func, dJxc)
end

function error!(pa::Param_error, grad=nothing)

	tgrid = pa.x.tgrid;
	acq = pa.x.acqgeom;
	fields = pa.x.fields;
	nss = acq.nss;

	# resample xr <-- x
	interp_spray!(pa.xr, pa.x, :interp)

	J = 0.0;
	for ifield=1:length(fields), iss=1:acq.nss
		yy=pa.y.d[iss, ifield]
		xx=pa.xr.d[iss,ifield]
		ww=pa.w.d[iss,ifield]

		dJxrr=pa.dJxr.d[iss, ifield]

		wav=pa.coup.ssf[iss,ifield]
		if(grad==:dJssf)
			dJwav=pa.dJssf[iss,ifield]
			dJwav[:]=0.0 # gradient is set to zero for all
		end
		
		for ir=1:acq.nr[iss]
			yyy=view(yy,:,ir)
			xxx=view(xx,:,ir)
			www=view(ww,:,ir)

			# apply source filter to xxx
			Conv.mod!(pa.paconvssf, :d, gf=xxx, wav=wav)
			xxxwav=pa.paconvssf.d

			if(grad==:dJx)
				# compute dJxc
				JJ = pa.func(pa.dJxc, xxxwav, yyy, www);

				dJxrrr=view(dJxrr,:,ir)
				# apply adjoint of source filter to dJxc
				Conv.mod!(pa.paconvssf, :gf, gf=dJxrrr, wav=wav, d=pa.dJxc)
				scale!(dJxrrr, inv(pa.ynorm))
			elseif(grad==:dJssf)
				# compute dJxc
				JJ = pa.func(pa.dJxc, xxxwav, yyy, www);

				# apply adjoint of source filter to dJxc
				Conv.mod!(pa.paconvssf, :wav, gf=xxx, d=pa.dJxc)
				for i in eachindex(dJwav)
					dJwav[i]+=pa.paconvssf.wav[i] # stack gradient of ssf over all receivers
				end
				scale!(dJwav, inv(pa.ynorm))
			else
				JJ = pa.func(nothing, xxxwav, yyy, www);
			end
			J += JJ 
		end
	end
	J /= pa.ynorm

	# spray dJx <-- dJxr
	interp_spray!(pa.dJx, pa.dJxr, :spray)

	(J == 0.0) && warn("misfit computed is zero")

	return J
end



end # module
