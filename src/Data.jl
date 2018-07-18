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
#using DeConv
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
Return a vec of data object sorted in the order
time, receivers, supersource, fields
"""
function Base.vec(data::TD)
	nr = data.acqgeom.nr;	nss = data.acqgeom.nss;	nt = data.tgrid.nx;
	v=Vector{Float64}()
	for ifield = 1:length(data.fields), iss = 1:nss
		dd=data.d[iss, ifield]
		append!(v,dd)
	end
	return v
end

"""
Fill with randn values
"""
function Base.randn!(data::TD)
	nr = data.acqgeom.nr;	nss = data.acqgeom.nss;	nt = data.tgrid.nx;
	for ifield = 1:length(data.fields), iss = 1:nss
		dd=data.d[iss, ifield]
		randn!(dd)
	end
	return data
end



"""
Opposite of method above, do not use in inner loops, allocations, not checked
"""
function Base.copy!(dataout::TD, v::Vector{Float64})
	nr = dataout.acqgeom.nr;	nss = dataout.acqgeom.nss;	nt = dataout.tgrid.nx;
	dout=getfield(dataout, :d)
	i0=0
	for iss=1:nss, ifield=1:length(dataout.fields)
		ddout=dout[iss,ifield]
		for ir = 1:nr[iss]
			for it in 1:nt
				ddout[it,ir]=v[i0+it]
			end
			i0+=nt
		end
	end
	return dataout
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
			# don't use explicit loop here, slow!
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


include("Data_misfits.jl")

#=
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

=#


end # module
