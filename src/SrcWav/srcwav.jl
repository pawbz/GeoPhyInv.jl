
using NamedArrays
using OrderedCollections

"""
Data type for the source related parameters during acquisiton.

# Fields
* `nss::Int64` : number of supersources
* `ns::Array{Int64}` : number of sources for each supersource
* `fields::Vector{Symbol}` : number of fields
* `wav::Array{Float64}` : wavelets in time domain
* `tgrid` : time grid 
"""
mutable struct SrcWavss
	tgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
	w::NamedArrays.NamedArray{Array{Float64,2},1,Array{Array{Float64,2},1},Tuple{OrderedCollections.OrderedDict{Symbol,Int64}}}
	ns::Int64
	#SrcWav(nss, ns, fields, wav, tgrid) = 
	#	any([
	#	  any([fields[ifield] ∉ [:P, :Vx, :Vz] for ifield in 1:length(fields)]),
	#	  length(fields) == 0,
	#	  broadcast(size,wav) ≠ [(length(tgrid),ns[iss]) for iss=1:nss, ifield=1:length(fields)]
	#	  ]) ? 
	#	error("error in SrcWav construction") : new(nss, ns, fields, wav, tgrid)
end

"""
Initialize an model with zeros, given mgrid and names.
"""
function Base.zero(::Type{SrcWavss}, tgrid, s::Srcs=Srcs(1), fields=[:P])
	nt=length(tgrid); nf=length(fields)
	return SrcWavss(tgrid,NamedArray([zeros(nt,s.n) for i in fields], (fields,)),
			s.n)
end
SrcWavss(tgrid, s::Srcs=Srcs(1), fields=[:P])=zero(SrcWavss, tgrid, s, fields)

"""
Print information about `SrcWav`
"""
function Base.show(io::Base.IO, src::SrcWavss)
	println(io,"\tSource")
#	println(io,"\t> number of supersources:\t",src.nss)
#	println(io,"\t> sources per supersource:\t","min\t",minimum(src.ns[:]), "\tmax\t", maximum(src.ns[:]))
	
#	freqmin, freqmax, freqpeak = freqs(src)
#	println(io,"\t> frequency:\t","min\t",freqmin, "\tmax\t",freqmax,"\tpeak\t",freqpeak)
#	println(io,"\t> time:\t","min\t",src.tgrid[1], "\tmax\t", src.tgrid[end])
#	println(io,"\t> samples:\t",length(src.tgrid))
end



#function Base.isapprox(src1::SrcWavss, src2::SrcWavss)
#	vec=([(src1.nss==src2.nss), (src1.fields==src2.fields), (src1.ns==src2.ns), 
#       		(isequal(src1.tgrid, src2.tgrid)),
#		])
#	return all(vec)
#end

function Base.copyto!(srco::SrcWavss, src::SrcWavss)
	@assert (srco.ns==src.n) && (srco.tgrid==src.tgrid)
	for i in 1:length(srco.w)
		copyto!(srco.w[i], src.w[i])
	end
end


function LinearAlgebra.rmul!(src::SrcWavss, x::Number)
	for i in 1:length(src.w)
		rmul!(src.w[i], x)
	end
end


SrcWav=Array{SrcWavss,1}

function Array{SrcWavss,1}(tgrid::StepRangeLen, ss::SSrcs=SSrcs(5), s::Vector{Srcs}=fill(Srcs(1),ss.n), fields=[:P])
	return [SrcWavss(tgrid,s[i], fields) for i in 1:ss.n]
end


#"""
#Allocate `SrcWav` with zeros depending on the acquisition geometry.
#"""
#function SrcWav_zeros(acqgeom::Geom,  fields::Vector{Symbol}, tgrid::StepRangeLen)
#	wavsrc = [zeros(length(tgrid),acqgeom.ns[iss]) for iss=1:acqgeom.nss, ifield=1:length(fields)] 
#	return SrcWav(acqgeom.nss, acqgeom.ns, fields, wavsrc, deepcopy(tgrid))
#end



"""
Return minimum, maximum and peak frequencies of `SrcWav`
"""
#function freqs(src::SrcWav)
#	freqmin = minimum([Utils.findfreq(src.wav[i,j][:,:],src.tgrid,attrib=:min) for i in 1:src.nss, j in 1:length(src.fields)])
#	freqmax = maximum([Utils.findfreq(src.wav[i,j][:,:],src.tgrid,attrib=:max) for i in 1:src.nss, j in 1:length(src.fields)])
#	freqpeak = mean([Utils.findfreq(src.wav[i,j][:,:],src.tgrid,attrib=:peak) for i in 1:src.nss, j in 1:length(src.fields)])
#	return freqmin, freqmax, freqpeak
#end


"""
Constructor for `SrcWav` data type.
Uses same source wavelet, i.e., `wav` for all sources and supersources

# Arguments

* `nss::Int64` : number of supersources
* `ns::Int64` : number of sources
* `fields::Vector{Symbol}` : number of fields the sources are exciting
* `wav::Array{Float64}` : a source wavelet that is used for all sources and supersources
* `tgrid` : time grid for the wavelet
"""
function update!(s::SrcWavss, fields::Vector{Symbol}, w::Array{Float64},)
	for f in fields
		for is in 1:s.ns
			ww=view(s.w[f],:,is)
			copyto!(ww,w)
		end
	end
	return s
end

function update!(s::SrcWav, fields::Vector{Symbol}, w::Array{Float64},)
	for ss in s
		update!(ss,fields,w)
	end
	return s
end




#=
"""
Constructor of `SrcWav`, which is typical for a input model such that 
the model has `nλ` wavelengths.

# Arguments

* `nss::Int64` : number of supersources
* `ns::Int64` : number of source per each supersource
* `fields::Vector{Symbol}` :

# Keyword Arguments

* `mod::Medium` :
* `nλ::Int64=10` : number of wavelengths in the mod
* `wav_func::Function=(fqdom, tgrid)->Utils.Wavelets.ricker(fqdom,tgrid)` : which wavelet to generate, see Utils.Wavelets.jl
* `tmaxfrac::Float64=1.0` : by default the maximum modelling time is computed using the average velocity and the diagonal distance of the model, use this fraction to increase of reduce the maximum time
"""
function update!(s::SrcWavss, 
		fields::Vector{Symbol},
		mod::Medium, 
		nλ::Int64,
		wav_func::Function=(fqdom, tgrid)->Utils.Wavelets.ricker(fqdom,tgrid),
		tmaxfrac::Float64=1.0
		)

	x=mod.mgrid[2]; z=mod.mgrid[1]
	# dominant wavelength using mod dimensions
	λdom=mean([(abs(x[end]-x[1])), (abs(z[end]-z[1]))])/real(nλ)
	# average P velocity
	vavg=mean(mod[:vp])

	fqdom = vavg/λdom

	# maximum distance the wave travels
	d = sqrt((x[1]-x[end]).^2+(z[1]-z[end]).^2)
	
	# use two-way maximum distance to get tmax
	tmax=2.0*d*inv(vavg)*tmaxfrac

	# choose sampling interval to obey max freq of source wavelet
	δmin = minimum([step(mod.mgrid[2]), step(mod.mgrid[1])])
	vpmax = mod.bounds[:vp][2]
	δt=0.5*δmin/vpmax

	# check if δt is reasonable
	#(δt > 0.1/fqdom) : error("decrease spatial sampling or nλ")

	tgrid=range(0.0, stop=tmax, step=δt)
	wav = wav_func(fqdom, tgrid);


	wavsrc = [hcat([wav_func(fqdom, tgrid) for is in 1:ns]...) for iss=1:nss, ifield=1:length(fields)] 
	src=SrcWav(nss, fill(ns, nss), fields, wavsrc, deepcopy(tgrid))
	print(src)
	# choose ricker waveletes of fdom
	return src
end
=#

#=

"""
Generate band-limited random source signals 
"""
function SrcWav_fixed_random(nss::Int64, ns::Int64, fields::Vector{Symbol}; 
			  distvec=[Normal() for iss in 1:nss],
			  sparsepvec=[1. for iss in 1:nss],
			  fmin::Float64=0.0, 
			  fmax::Float64=0.0, 
			  tgrid::StepRangeLen=nothing, 
			  tmaxfrac::Float64=1.0)
	wavsrc = [repeat(zeros(length(tgrid)),inner=(1,ns)) for iss=1:nss, ifield=1:length(fields)] 
	for ifield in 1:length(fields), iss in 1:nss, is in 1:ns
		wavsrc[iss, ifield][:,is] = Utils.get_tapered_random_tmax_signal(tgrid, fmin=fmin, fmax=fmax, tmaxfrac=tmaxfrac, dist=distvec[iss], 
								 sparsep=sparsepvec[iss])
	end
	src=SrcWav(nss, fill(ns, nss), fields, wavsrc, deepcopy(tgrid))
	print(src)
	return src
end


"""
Function that returns SrcWav after time reversal
"""
function SrcWav_tr(src::SrcWav)
	return SrcWav(src.nss,src.ns,src.fields,
	    [reverse(src.wav[i,j],dims=1) for i in 1:src.nss, j in 1:length(src.fields)],deepcopy(src.tgrid))
end



"""
Pad `SrcWav` 
tgrids should be same in all SrcWav
"""
function SrcWav_uspos(src::Vector{SrcWav}, acq::Vector{Geom})
	np = length(src) == length(acq) ? length(src) : error("unequal sizez")

	# unique source positions
	nus=Geom_get(acq,:nus) 
	uspos=Geom_get(acq,:uspos)
	# all zeros for all unique positions
	wavout = [[zeros(src[ip].length(tgrid),nus) for iss=1:src[ip].nss, ifield=1:length(src[ip].fields)] for ip=1:np]
	# fill source wavelets when necessary
	for ip=1:np
		for ifield=1:length(src[ip].fields), iss=1:src[ip].nss, is=1:acq[ip].ns[iss]
			is0=findall([[uspos[1][i]-acq[ip].sz[iss][is],
		      uspos[2][i]-acq[ip].sx[iss][is]] == [0., 0.,] for i in 1:nus])

			wavout[ip][iss, ifield][:,is0] = src[ip].wav[iss, ifield][:,is] 
		end
	end
	# output src
	return [SrcWav(src[ip].nss,fill(nus, src[ip].nss),
	     src[ip].fields,wavout[ip],src[ip].tgrid) for ip=1:np]
end




"""
return a vector of the order 

"""
function SrcWav_getvec(src::Vector{SrcWav}, field::Symbol)
	np = length(src);
	vect = [getfield(src[ip],field) for ip=1:np]
	return vec(hcat(hcat(vect...)...))
end

=#
