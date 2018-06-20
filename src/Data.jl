__precompile__()

"""
This module defines the data types related to seismic data:
* `TD` : time domain representation

It also provides methods that apply source and receiver filters onto 
seismic data.
"""
module Data

using Grid
using Conv
using Interpolation
using DeConv
using Misfits
using Signals
import JuMIT.Acquisition
import JuMIT.Coupling
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
	interp_spray!(data, dataout)
	return dataout
end


"""
Method to resample data in time.
Can reduce allocations =========

# Arguments

* `data::TD` : input data of type `TD`
* `dataout::TD` : preallocated data of type `TD` that is modified
"""
function interp_spray!(data::TD, dataout::TD, attrib=:interp, Battrib=:B1; pa=nothing)
	# check if datasets are similar
	nss = data.acqgeom.nss
	nr = data.acqgeom.nr
	xin=data.tgrid.x
	xout=dataout.tgrid.x
	if(pa===nothing)
		pa=Interpolation.Param([xin], [xout], :B1)
	end
	for ifield = 1:length(data.fields), iss = 1:nss
		dat=data.d[iss,ifield]
		dato=dataout.d[iss,ifield]
		for ir = 1:nr[iss]
			din=view(dat,:,ir)
			dout=view(dato,:,ir)
	 		Interpolation.interp_spray!(din, dout, pa, attrib)
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
Returns the variance of data

"""
function Base.var(data1::TD)
	σ=0.0
	μ=mean(data1)
	n=0
	for ifield = 1:length(data1.fields), iss = 1:data1.acqgeom.nss 
		for ir = 1:data1.acqgeom.nr[iss], it = 1:data1.tgrid.nx
			n += 1
			σ += (data1.d[iss, ifield][it, ir]-μ)^2 
		end
	end
	return σ*inv(n)
end

function Base.mean(data1::TD)
	n=0
	μ=0.0
	for ifield = 1:length(data1.fields), iss = 1:data1.acqgeom.nss 
		for ir = 1:data1.acqgeom.nr[iss], it = 1:data1.tgrid.nx
			n += 1
			μ += data1.d[iss, ifield][it, ir]
		end
	end
	return μ*inv(n)
end

function addnoise!(dataN::TD, data::TD, SNR)

	σx=var(data)

	σxN=sqrt(σx^2*inv(10^(SNR/10.)))
	
	# factor to be multiplied to each scalar
	α=sqrt(σxN)
	for ifield = 1:length(data.fields), iss = 1:data.acqgeom.nss 
		for ir = 1:data.acqgeom.nr[iss], it = 1:data.tgrid.nx

			dataN.d[iss, ifield][it, ir] = dataN.d[iss, ifield][it, ir] + α*randn()
		end
	end
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
	#	Signals.DSP.fast_filt!(s.d[iss, ifield][:, ir],r.d[iss, ifield][:, ir],
#		 w.rf[iss, ifield][:,ir],:s)
		# source coupling
		sv=s.d[iss, ifield][:, ir]
		rv=r.d[iss, ifield][:, ir]
		wv=w.ssf[iss, ifield]
		Signals.DSP.fast_filt!(sv, rv, wv, attrib)
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
	twin[itlim[1] : itlim[2]] = Signals.DSP.taper(ones(itlim[2]-itlim[1]+1),ttaperperc) 

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
	xrc::TD # modelled data after resampling and convolution with source coupling
	dJxr::TD # gradient w.r.t. xr
	dJxrc::TD # gradient w.r.t. xr
	dJx::TD # gradient w.r.t. x
	dJssf::Matrix{Vector{Float64}} # gradient w.r.t source filters
	ynorm::Float64 # normalize functional with energy of y
	coup::Coupling.TD
	paconvssf::Conv.P_conv{Float64,1} # convolutional model for source filters
	paconvrf::Conv.P_conv{Float64,1} # convolutional model for receiver filters
	pacse::Matrix{Conv.P_misfit_xcorr}
	func::Matrix{Function} # function to compute misfit for every ss and fields
	painterp::Interpolation.Param{Float64}
end

function Param_error(x,y;w=nothing, coup=nothing, func_attrib=:cls)

	if(coup===nothing)
		 coup=Coupling.TD_delta(y.tgrid,[0.1,0.1],[0.0, 0.0], [:P], y.acqgeom)
	end

	paconvssf=Conv.P_conv(ssize=[coup.tgridssf.nx], 
		      dsize=[y.tgrid.nx], 
		      gsize=[y.tgrid.nx], 
		      slags=coup.ssflags, 
		      dlags=[y.tgrid.nx-1, 0], 
		      glags=[y.tgrid.nx-1, 0])

	paconvrf=Conv.P_conv(ssize=[coup.tgridrf.nx], 
		      dsize=[y.tgrid.nx], 
		      gsize=[y.tgrid.nx], 
		      slags=coup.rflags, 
		      dlags=[y.tgrid.nx-1, 0], 
		      glags=[y.tgrid.nx-1, 0])

	dJssf=deepcopy(coup.ssf)

	painterp=Interpolation.Param([x.tgrid.x], [y.tgrid.x], :B1)

	if(w===nothing) 
		w=Data.TD_ones(y.fields,y.tgrid,y.acqgeom)
	end
	!(isapprox(w,y)) && error("weights have to be similar to y")

	!(isequal(x.acqgeom, y.acqgeom)) && error("observed and modelled data should have same acqgeom")
	!(isequal(x.fields, y.fields)) && error("observed and modelled data should have same fields")

	xr=TD_zeros(y)
	xrc=TD_zeros(y)
	dJxr=TD_zeros(y)
	dJxrc=TD_zeros(y)
	dJx=TD_zeros(x)

	ynorm = dot(y, y)
	(ynorm == 0.0) && error("y cannot be zero")

	dJxc=zeros(y.tgrid.nx)

	if(func_attrib==:cls)
		pacse=[Conv.P_misfit_xcorr(1, 1,y=zeros(1,1)) for i in 1:2, j=1:2] # dummy
		func=[(dJx,x)->Misfits.error_squared_euclidean!(dJx,x,y.d[iss,ifield],w.d[iss,ifield]) for iss in 1:y.acqgeom.nss, ifield=1:length(y.fields)]
	elseif(func_attrib==:xcorrcls)
		pacse=[Conv.P_misfit_xcorr(y.tgrid.nx, y.acqgeom.nr[iss],y=y.d[iss,ifield]) for iss in 1:y.acqgeom.nss, ifield=1:length(y.fields)]
		func=[(dJx,x)->Misfits.error_corr_squared_euclidean!(dJx,x,pacse[iss,ifield]) for iss in 1:y.acqgeom.nss, ifield=1:length(y.fields)]
	end

	pa=Param_error(x,y,w,xr,xrc,dJxr,dJxrc,
		    dJx,dJssf,ynorm,coup,paconvssf, paconvrf, pacse, func, painterp)


	return pa
end

function error!(pa::Param_error, grad=nothing)

	tgrid = pa.x.tgrid;
	acq = pa.x.acqgeom;
	fields = pa.x.fields;
	nss = acq.nss;

	# resample xr <-- x
	interp_spray!(pa.x, pa.xr, :interp, pa=pa.painterp)

	xr=pa.xr
	xrc=pa.xrc
	y=pa.y


	J = 0.0;
	for ifield=1:length(fields), iss=1:acq.nss
		xrr=xr.d[iss,ifield]
		xrcc=xrc.d[iss,ifield]
		wav=pa.coup.ssf[iss,ifield]

		nt=size(xrr,1)
		nr=size(xrr,2)

		copy!(pa.paconvssf.s,wav)
		# xrc <- xr apply source filter to xr
		for ir=1:nr
			xrrr=view(xrr,:,ir)
			for i in 1:nt
				pa.paconvssf.d[i]=xrr[i,ir]
			end
			Conv.mod!(pa.paconvssf, :d, g=xrrr, s=wav)
			for i in 1:nt
				xrcc[i,ir]=pa.paconvssf.d[i]
			end
		end

		if(grad===nothing)
			JJ=pa.func[iss,ifield](nothing,  xrcc)
		elseif((grad==:dJx) || (grad==:dJssf))
			dJxrcc=pa.dJxrc.d[iss,ifield]
			JJ=pa.func[iss,ifield](dJxrcc,  xrcc)
		else
			error("invalid grad")
		end

		# dJxr <- dJxrc  apply adjoint of source filter to dJxc
		if(grad==:dJx)
			dJxrr=pa.dJxr.d[iss,ifield]
			dJxrcc=pa.dJxrc.d[iss,ifield]
			copy!(pa.paconvssf.s,wav)
			for ir=1:nr
				for i in 1:nt
					pa.paconvssf.d[i]=dJxrcc[i,ir]
				end
				Conv.mod!(pa.paconvssf, :g)
				for i in 1:nt
					dJxrr[i,ir]=pa.paconvssf.g[i]
				end
			end
			scale!(dJxrr, inv(pa.ynorm))  # take care of scale later in the functional
		end

		# dJssf <- dJxrc 
		if(grad==:dJssf)
			dJwav=pa.dJssf[iss,ifield]
			dJwav[:]=0.0 # gradient is set to zero for all
			for ir=1:acq.nr[iss]
				xrrr=view(xrr,:,ir)
				dJxrccc=view(dJxrcc,:,ir)
				Conv.mod!(pa.paconvssf, :s, g=xrrr, d=dJxrccc)
				for i in eachindex(dJwav)
					dJwav[i]+=pa.paconvssf.s[i] # stack gradient of ssf over all receivers
				end
			end
			scale!(dJwav, inv(pa.ynorm))
		end
		J += JJ 
	end
	J /= pa.ynorm

	# spray dJx <-- dJxr
	interp_spray!(pa.dJx, pa.dJxr, :spray, pa=pa.painterp)

	(J == 0.0) && warn("misfit computed is zero")

	return J

	J = 0.0;
	for ifield=1:length(fields), iss=1:acq.nss
		pacsee=pa.pacse[iss,ifield]
		error_corr_squared_euclidean!(dfdx,  x; pa=pacsee)

	end
end


function DDeConv(d::TD, wav::AbstractVector{Float64}, ϵ=1e-2)

	dout=TD_zeros(d)
	ntd=dout.tgrid.nx
	wavv=deepcopy(wav);

	paD=DeConv.ParamD(ntd=ntd,nts=length(wav), s=wavv)

	paD.ϵ=ϵ

	DDeConv!(dout, d, paD)
	return dout
end


function DDeConv!(dataout::TD, data::TD, paD)
	nr = data.acqgeom.nr;	nss = data.acqgeom.nss;	nt = data.tgrid.nx;
	for ifield = 1:length(data.fields), iss = 1:nss
		dd=data.d[iss, ifield]
		ddo=dataout.d[iss, ifield]
		for ir = 1:nr[iss]
			for it in 1:nt
				paD.d[it]=dd[it,ir]
			end
			DeConv.mod!(paD)
			for it in 1:nt
				ddo[it,ir]=paD.g[it]
			end
		end
	end
	return dataout
end



end # module
