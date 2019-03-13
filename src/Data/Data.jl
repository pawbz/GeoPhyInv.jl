"""
This module defines the data types related to seismic data:
* `TD` : time domain representation

It also provides methods that apply source and receiver filters onto 
seismic data.
"""
module Data

using Conv
using Interpolation
#using DeConv
using Misfits
using Signals
using LinearAlgebra
using Statistics
import JuMIT.Acquisition
import JuMIT.Coupling
using DSP
using Random

"""
Time domain representation of Seismic Data.

# Fields

* `d::Array{Array{Float64,2},2}` : data 
* `fields::Vector{Symbol}` : components recorded at each receiver
* `tgrid` : grid to represent time
* `acqgeom::Acquisition.Geom` : acquisition geometry used to generate the data
"""
mutable struct TD
	d::Matrix{Matrix{Float64}}
	fields::Vector{Symbol}
	tgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
	acqgeom::Acquisition.Geom
	"adding conditions that are to be false while construction"
	TD(d, fields, tgrid, acqgeom) = 
		any([
		  any([fields[ifield] ∉ [:P, :Vx, :Vz] for ifield in 1:length(fields)]),
		  length(fields) == 0,
		  broadcast(size,d) != [(length(tgrid),acqgeom.nr[iss]) for iss=1:acqgeom.nss, ifield=1:length(fields)]
		  ]) ? 
		error("error in TD construction") : new(d, fields, tgrid, acqgeom)

end


"""
Method to resample data in time.

# Arguments

* `data::TD` : input data of type `TD`
* `tgrid` : resampling in time according to this time grid

# Return

* data after resampling as `TD`
"""
function interp(data::TD,
		tgrid::StepRangeLen
		)
	nss = data.acqgeom.nss
	nr = data.acqgeom.nr
	dataout = TD(
	      [zeros(length(tgrid),data.acqgeom.nr[iss]) for iss=1:nss, ifield=1:length(data.fields)],
	      data.fields,tgrid,data.acqgeom)
	interp_spray!(data, dataout)
	return dataout
end

function taper!(data::TD, perc=0.0; bperc=perc,eperc=perc)
	nr = data.acqgeom.nr;	nss = data.acqgeom.nss;	nt = length(data.tgrid);
	for ifield = 1:length(data.fields), iss = 1:nss
		dd=data.d[iss, ifield]
		Signals.DSP.taper!(dd,bperc=bperc,eperc=eperc)
	end
	return data
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
	xin=data.tgrid
	xout=dataout.tgrid
	if(pa===nothing)
		pa=Interpolation.Kernel([xin], [xout], :B1)
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
* `tgrid` : time domain grid
* `acqgeom::Acquisition.Geom` : acquisition geometry

# Return

* data with zeros as `TD`
"""
function TD_zeros(fields::Vector{Symbol}, tgrid::StepRangeLen, acqgeom::Acquisition.Geom)
	return TD([zeros(length(tgrid),acqgeom.nr[iss]) for iss=1:acqgeom.nss, ifield=1:length(fields)],fields,
	   deepcopy(tgrid),deepcopy(acqgeom)) 
end
function Base.fill!(data::TD, k::Float64)
	for iss=1:data.acqgeom.nss, ifield=1:length(data.fields)
		fill!(data.d[iss,ifield],k) 
	end
end
function TD_zeros(d::TD)
	return TD([zeros(length(d.tgrid),d.acqgeom.nr[iss]) for iss=1:d.acqgeom.nss, ifield=1:length(d.fields)],d.fields,
	   deepcopy(d.tgrid),deepcopy(d.acqgeom)) 
end


"Same as `TD_zeros`, except for returning ones"
function TD_ones(fields::Vector{Symbol}, tgrid::StepRangeLen, acqgeom::Acquisition.Geom) 
	return TD([ones(length(tgrid),acqgeom.nr[iss]) for iss=1:acqgeom.nss, ifield=1:length(fields)],
	   fields,deepcopy(tgrid),deepcopy(acqgeom)) 
end



"""
Time reverse the records of each receiver in `TD` 

# Arguments

* `data::TD` : input data that is modified
"""
function TD_tr!(data::TD)
	data.d = copy([reverse(data.d[i,j],dims=1) for i in 1:data.acqgeom.nss, j in 1:length(data.fields)]);
end


function addnoise!(dataN::TD, data::TD, SNR)

	σx=Statistics.var(data)

	σxN=sqrt(σx^2*inv(10^(SNR/10.)))
	
	# factor to be multiplied to each scalar
	α=sqrt(σxN)
	for ifield = 1:length(data.fields), iss = 1:data.acqgeom.nss 
		for ir = 1:data.acqgeom.nr[iss], it = 1:length(data.tgrid)

			dataN.d[iss, ifield][it, ir] = dataN.d[iss, ifield][it, ir] + α*Random.randn()
		end
	end
end



# include rest of the files
for file in ["base", "statistics", "processing", "misfits", "weights"]   
	fn=joinpath(@__DIR__, string(file,".jl"))
	include(fn)
	#Revise.track(@__MODULE__,fn)
end


#=
function DDeConv(d::TD, wav::AbstractVector{Float64}, ϵ=1e-2)

	dout=TD_zeros(d)
	ntd=length(dout.tgrid)
	wavv=deepcopy(wav);

	paD=DeConv.ParamD(ntd=ntd,nts=length(wav), s=wavv)

	paD.ϵ=ϵ

	DDeConv!(dout, d, paD)
	return dout
end


function DDeConv!(dataout::TD, data::TD, paD)
	nr = data.acqgeom.nr;	nss = data.acqgeom.nss;	nt = length(data.tgrid);
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
