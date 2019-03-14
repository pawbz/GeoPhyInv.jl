module Models


using DataFrames
using CSV
using Statistics
using LinearAlgebra
using Distributions
using Random
#using ImageFiltering
import JuMIT.Smooth
import JuMIT.Interpolation

"""
Store reference model parameters
"""
mutable struct Seismic_ref{T<:Real}
	vp::T
	vs::T
	ρ::T
	K::T
	μ::T
	KI::T
	μI::T
	ρI::T
end


"""
Data type fo represent a seismic model.
A contrast function for a model m is given by ``χ(m) = \frac{m}{m0}-1``.

# Fields

* `vp0::Vector{Float64}` : [vpmin, vpmax]
* `vs0::Vector{Float64}` : [vsmin, vsmax]
* `ρ0::Vector{Float64}` : [ρmin, ρmax]
* `χvp::Array{Float64,2}` : two-dimensional contrast model (χ) for vp, for e.g., zeros(length(mgrid[1]), length(mgrid[2]))
* `χvs::Array{Float64}` : two-dimensional contrast model (χ) for vs, for e.g., zeros(length(mgrid[1]), length(mgrid[2]))
* `χρ::Array{Float64}` : two-dimensional contrast model (χ) for density, for e.g., zeros(length(mgrid[1]), length(mgrid[2]))
* `mgrid` : array of ranges to determine the dimensions of models
"""
mutable struct Seismic
	vp0::Vector{Float64}
	vs0::Vector{Float64}
	ρ0::Vector{Float64}
	χvp::Array{Float64,2}
	χvs::Array{Float64,2}
	χρ::Array{Float64,2}
	mgrid::Vector{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}}
	# all derived fields from previous
	K0::Vector{Float64}
	μ0::Vector{Float64}
	KI0::Vector{Float64}
	μI0::Vector{Float64}
	ρI0::Vector{Float64}
	ref::Seismic_ref{Float64}
	"adding conditions that are to be false while construction"
	Seismic(vp0,vs0,ρ0,χvp,χvs,χρ,mgrid,K0,μ0,KI0,μI0,ρI0,ref) = 
		any([
		     any(vp0.<0.0), 
		     any(vs0.<0.0),
		     any(ρ0.<0.0),
		     #all([all(vp0 .≠ 0.0), any(χ(χvp,vp0,-1) .< vp0[1])]), # check vp bounds
		     #all([all(vp0 .≠ 0.0), any(χ(χvp,vp0,-1) .> vp0[2])]), # check vp bounds
		     #all([all(vs0 .≠ 0.0), any(χ(χvs,vs0,-1) .< vs0[1])]), # check vs bounds
		     #all([all(vs0 .≠ 0.0), any(χ(χvs,vs0,-1) .> vs0[2])]), # check vs bounds
		     #all([all(ρ0 .≠ 0.0), any(χ(χρ,ρ0,-1) .< ρ0[1])]), # check ρ bounds
		     #all([all(ρ0 .≠ 0.0), any(χ(χρ,ρ0,-1) .> ρ0[2])]), # check ρ bounds
		     size(χvp) != (length(mgrid[1]), length(mgrid[2])), # dimension check
		     size(χvs) != (length(mgrid[1]), length(mgrid[2])), # dimension check
		     size(χρ) != (length(mgrid[1]), length(mgrid[2])) # dimension check
		    ]) ? 
		       error("error in Seismic construction") : new(vp0,vs0,ρ0,χvp,χvs,χρ,mgrid,K0,μ0,KI0,μI0,ρI0,ref)
end

# check whether a seismic model is bounded
function isbounded(mod::Seismic)
	@warn "need to be implemented"
	result=true
end

# method to create Seismic object without using derived fields
function Seismic(vp0, vs0, ρ0, χvp, χvs, χρ, 
		 mgrid::Vector{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}})
	mod=Seismic(vp0, vs0, ρ0, χvp, χvs, χρ, mgrid,
	    zeros(2), zeros(2), zeros(2), zeros(2), zeros(2),
	    Seismic_ref(zeros(8)...))
	update_derived!(mod)
	return mod
end

# update the derived fields of the Seisimc model
function update_derived!(mod::Seismic)
	"note that: mean(K0) ≠ 1/mean(KI0)"
	mod.K0 = mod.vp0 .* mod.vp0 .* mod.ρ0; 
	mod.KI0 = reverse(inv.(mod.K0),dims=1);
	ρ0 = mod.ρ0; 
	mod.ρI0 = reverse(inv.(ρ0),dims=1);
	mod.μ0 = mod.vs0 .* mod.vs0 .* mod.ρ0;
	mod.μI0 = reverse(inv.(mod.μ0),dims=1);
	mod.ref = Seismic_ref(Statistics.mean.([mod.vp0,mod.vs0,mod.ρ0,mod.K0,mod.μ0,mod.KI0,mod.μI0,mod.ρI0])...)
end

"""
Print information about `Seismic`
"""
function Base.print(mod::Seismic, name::String="")
	println("\tSeismic Model:\t",name)
	println("\t> number of samples:\t","x\t",length(mod.mgrid[2]),"\tz\t",length(mod.mgrid[1]))
	println("\t> sampling intervals:\t","x\t",step(mod.mgrid[2]),"\tz\t",step(mod.mgrid[1]))
	println("\t> vp:\t","min\t",minimum(Seismic_get(mod,:vp)),"\tmax\t",maximum(Seismic_get(mod,:vp)))
	println("\t> vp bounds:\t","min\t",mod.vp0[1],"\tmax\t",mod.vp0[2])
	println("\t> ρ:\t","min\t",minimum(Seismic_get(mod,:ρ)),"\tmax\t",maximum(Seismic_get(mod,:ρ)))
	println("\t> ρ bounds:\t","min\t",mod.ρ0[1],"\tmax\t",mod.ρ0[2])
end

"""
Return medium property bounds based on maximum and minimum values of the array and frac.
The bounds cannot be less than zero
"""
function bounds(mod::Array{Float64}, frac::Float64=0.1)
	any(mod .< 0.0) && error("model values less than zero")
	bounds=zeros(2)
	bound=frac*Statistics.mean(mod)
	bounds[1] = ((minimum(mod) - bound)<0.0) ? 0.0 : (minimum(mod) - bound)
	bounds[2] = maximum(mod)+bound
	return bounds
end


function adjust_bounds!(mod::Seismic,mod0::Seismic)
	adjust_bounds!(mod,mod0.vp0,mod0.vs0,mod0.ρ0)
end

function adjust_bounds!(mod::Seismic,vp0::Vector{T},vs0::Vector{T},ρ0::Vector{T}) where {T<:Real}
	copyto!(mod.vp0, Float64.(vp0))
	copyto!(mod.vs0, Float64.(vs0))
	copyto!(mod.ρ0, Float64.(ρ0))
	update_derived!(mod)
end

"""
Adjust the bounds and hence the reference values.
Since the reference values are adjust the χ fields should also be changed
"""
function adjust_bounds!(mod::Seismic, frac::Float64=0.1)
	for f in [:χvp, :χρ, :χvs]
		f0 = Symbol((replace("$(f)", "χ" => "")),0)
		fm = Symbol((replace("$(f)", "χ" => "")))
		m = Seismic_get(mod, fm)
		m0 = bounds(m, frac)
		setfield!(mod, f0, m0)
		setfield!(mod, f, χ(m, mean(m0), 1))
	end
	update_derived!(mod)
	return mod
end

"""
Return `Seismic` with zeros everywhere;
this method is used for preallocation.

# Arguments
* `mgrid` : used for sizes of χ fields 
"""
function Seismic_zeros(mgrid::Vector{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}})
	return Seismic(fill(0.0,2), fill(0.0,2), fill(0.0,2), 
		zeros(length(mgrid[1]), length(mgrid[2])), zeros(length(mgrid[1]), length(mgrid[2])),
		zeros(length(mgrid[1]), length(mgrid[2])), deepcopy(mgrid))
end

function Base.fill!(mod::Seismic, k::Float64=0.0)
	fill!(mod.χvp,k)
	fill!(mod.χρ,k)
	fill!(mod.χvs,k)
	return mod
end

"""
Return true if a `Seismic` model is just allocated with zeros.
"""
function Base.iszero(mod::Seismic)
	return any((mod.vp0 .* mod.ρ0) .== 0.0) ? true : false
end

"""
Return a similar model to the input model, used for allocation.
"""
function Base.similar(mod::Seismic)
	modout=deepcopy(mod)
	fill!(modout, 0.0)
	return modout
end


"Compare if two `Seismic` models are equal"
function Base.isequal(mod1::T, mod2::T) where {T<:Union{Seismic,Seismic_ref}}
	fnames = fieldnames(T)
	vec=[(isequal(getfield(mod1, name),getfield(mod2, name))) for name in fnames]
	return all(vec)
end

"""
Return if two `Seismic` models have same dimensions and bounds.
"""
function Base.isapprox(mod1::Seismic, mod2::Seismic)
	vec=([(mod1.vp0==mod2.vp0), (mod1.vs0==mod2.vs0), (mod1.ρ0==mod2.ρ0), 
       		(size(mod1.χvp)==size(mod2.χvp)), 
       		(size(mod1.χvs)==size(mod2.χvs)), 
       		(size(mod1.χρ)==size(mod2.χρ)), 
		isequal(mod1.mgrid, mod2.mgrid),
		])
	return all(vec)
end

"""
Copy for `Seismic` models. The models should have same bounds and sizes.
Doesn't allocate any memory.
"""
function Base.copyto!(modo::Seismic, mod::Seismic)
	if(isapprox(modo, mod))
		for f in [:χvp, :χρ, :χvs]
			modof=getfield(modo,f)
			modf=getfield(mod,f)
			copyto!(modof,modf)
		end
		return modo
	else
		error("attempt to copy dissimilar models")
	end
end

function Seismic_get(mod::Seismic, attrib::Symbol)
	# allocate
	rout=zeros(length(mod.mgrid[1]), length(mod.mgrid[2]))
	Seismic_get!(rout, mod, [attrib])
	return rout
end

"""
Get other dependent model parameters of a seismic model
that are not present in `Seismic`.

* `:ρI` : inverse of density
* `:Zp` : P-wave impedance
"""
function Seismic_get!(x, mod::Seismic, attribvec::Vector{Symbol})
	ρ=mod.χρ; χ!(ρ, mod.ref.ρ,-1) # undo it later
	vp=mod.χvp; χ!(vp, mod.ref.vp,-1) # undo it later
	vs=mod.χvs; χ!(vs, mod.ref.vs,-1) # undo it later
	nznx=length(ρ)
	(length(x)≠(count(attribvec.≠ :null)*nznx)) &&  error("size x")
	i0=0; 
	for attrib in attribvec
		if(attrib ≠:null)
		if(attrib == :ρI)
			@inbounds for i in 1:nznx; x[i0+i]=inv(ρ[i]); end
		elseif(attrib == :ρ)
			@inbounds for i in 1:nznx; x[i0+i]=ρ[i]; end
		elseif(attrib == :vp)
			@inbounds for i in 1:nznx; x[i0+i]=vp[i]; end
		elseif(attrib == :vs)
			@inbounds for i in 1:nznx; x[i0+i]=vs[i]; end
		elseif(attrib == :χρI)
			@inbounds for i in 1:nznx; x[i0+i]=χ(inv(ρ[i]),mod.ref.ρI,1); end
		elseif(attrib == :χρ)
			@inbounds for i in 1:nznx; x[i0+i]=χ(ρ[i],mod.ref.ρ,1); end
		elseif(attrib == :χvp)
			@inbounds for i in 1:nznx; x[i0+i]=χ(vp[i],mod.ref.vp,1); end
		elseif(attrib == :χK)
			@inbounds for i in 1:nznx; x[i0+i]=χ(vp[i]*vp[i]*ρ[i],mod.ref.K,1); end
		elseif(attrib == :χμ)
			@inbounds for i in 1:nznx; x[i0+i]=χ(vs[i]*vs[i]*ρ[i],mod.ref.μ,1); end
		elseif(attrib == :χKI)
			@inbounds for i in 1:nznx; x[i0+i]=χ(inv(vp[i]*vp[i]*ρ[i]),mod.ref.KI,1); end
		elseif(attrib == :KI)
			@inbounds for i in 1:nznx; x[i0+i]=inv(vp[i]*vp[i]*ρ[i]); end
		elseif(attrib == :K)
			@inbounds for i in 1:nznx; x[i0+i]=vp[i]*vp[i]*ρ[i]; end
		elseif(attrib == :Zp)
			@inbounds for i in 1:nznx; x[i0+i]=vp[i]*ρ[i]; end
		else
			error("invalid attrib")
		end
		i0+=nznx
		end
	end
	χ!(ρ, mod.ref.ρ,1) 
	χ!(vp, mod.ref.vp,1)
	χ!(vs, mod.ref.vs,1)
	return x
end

include("reparameterize.jl")

include("chainrule.jl")

"""
Return dimensionless contrast model parameter using the
reference value.

# Arguments
* `mod::Array{Float64}` : subsurface parameter
* `mod0::Vector{Float64}` : reference value is mean of this vector
* `flag::Int64=1` : 
"""
function χ(mod::AbstractArray{T}, mod0::T, flag::Int64=1) where {T<:Real}
	mod_out=copy(mod)
	χ!(mod_out, mod0, flag)
	return mod_out
end
function χ!(mod::AbstractArray{T}, mod0::T, flag::Int64=1) where {T<:Real}
	m0 = mod0
	if(flag == 1)
		if(iszero(m0))
			return nothing
		else
			@inbounds for i in eachindex(mod)
				mod[i] = ((mod[i] - m0) * inv(m0))
			end
			return mod 
		end
	elseif(flag == -1)
		@inbounds for i in eachindex(mod)
			mod[i] = mod[i] * m0 + m0
		end
		return mod 
	end
end

function χ(m::T, m0::T, flag::Int64=1) where {T<:Real}
	if(flag==1)
		return (m-m0) * inv(m0)
	else(flag==-1)
		return (m*m0+m0)
	end
end

"""
Gradients
Return contrast model parameter using the
reference value.
"""
function χg(mod::AbstractArray{T}, mod0::T, flag::Int64=1) where {T<:Real}
	mod_out=copy(mod)
	χg!(mod_out, mod0, flag)
	return mod_out
end
function χg!(mod::AbstractArray{T}, mod0::T, flag::Int64=1) where {T<:Real}
	if(flag == 1)
		rmul!(mod, mod0)
	elseif(flag == -1)
		m0 = inv(mod0)
		rmul!(mod, m0)
	end
	return mod
end

function χg(m::T, m0::T, flag::Int64=1) where {T<:Real}
	if(flag==1)
		return m*m0
	else(flag==-1)
		return m*inv(m0)
	end
end




"""
Add features to a model.

# Arguments
* `mod::Seismic` : model that is modified

# Keyword Arguments

* `point_loc::Vector{Float64}=[0., 0.,]` : approx location of point pert.
* `point_pert::Float64=0.0` : perturbation at the point scatterer
* `ellip_loc::Vector{Float64}=nothing` : location of center of perturbation, [z, x]
* `ellip_rad::Float64=0.0` : radius of circular perturbation
* `ellip_pert::Float64=0.1` : perturbation inside a circle
* `rect_loc::Array{Float64}=nothing` : rectangle location, [zmin, xmin, zmax, xmax]
* `rect_pert::Float64=0.1` : perturbation in a rectangle
* `randn_pert::Float64=0.0` : percentage of reference values for additive random noise
* `fields::Vector{Symbol}=[:χvp,:χρ,:χvs]` : which fields are to be modified?
* `onlyin` : `mod` is modified only when field values are in these ranges 
"""
function Seismic_addon!(mod::Seismic; 
		       point_loc=[0., 0.,],
		       point_pert::Float64=0.0,
		       ellip_loc=[0., 0.,],
		       ellip_rad=0.0,
		       ellip_pert::Float64=0.0,
		       ellip_α=0.0,
		       rect_loc=[0., 0., 0., 0.,],
		       rect_pert::Float64=0.0,
		       constant_pert::Float64=0.0,
		       randn_perc::Real=0.0,
		       fields::Vector{Symbol}=[:χvp,:χρ,:χvs],
		       onlyin::Vector{Vector{Float64}}=[[typemin(Float64), typemax(Float64)] for i in fields]
		       )
	rect_loc=convert.(Float64,rect_loc);
	ellip_loc=convert.(Float64,ellip_loc);

	temp = zeros(length(mod.mgrid[1]), length(mod.mgrid[2]))

	ipointlocx = Interpolation.indminn(mod.mgrid[2], Float64(point_loc[2]), 1)[1] 
	ipointlocz = Interpolation.indminn(mod.mgrid[1], Float64(point_loc[1]), 1)[1] 
	temp[ipointlocz, ipointlocx] += point_pert

	if(!(ellip_pert == 0.0))
		α=ellip_α*pi/180.
		# circle or ellipse
		rads= (length(ellip_rad)==1) ? [ellip_rad[1],ellip_rad[1]] : [ellip_rad[1],ellip_rad[2]]

		temp += [(((((mod.mgrid[2][ix]-ellip_loc[2])*cos(α)+(mod.mgrid[1][iz]-ellip_loc[1])*sin(α))^2*inv(rads[1]^2) + 
	      ((-mod.mgrid[1][iz]+ellip_loc[1])*cos(α)+(mod.mgrid[2][ix]-ellip_loc[2])*sin(α))^2*inv(rads[2]^2)) <= 1.) ? ellip_pert : 0.0)  for 
	   		iz in 1:length(mod.mgrid[1]), ix in 1:length(mod.mgrid[2]) ]
	end
	if(!(rect_pert == 0.0))
		temp += [
			(((mod.mgrid[2][ix]-rect_loc[4]) * (mod.mgrid[2][ix]-rect_loc[2]) < 0.0) & 
			((mod.mgrid[1][iz]-rect_loc[3]) * (mod.mgrid[1][iz]-rect_loc[1]) < 0.0)) ?
			rect_pert : 0.0  for
			iz in 1:length(mod.mgrid[1]), ix in 1:length(mod.mgrid[2]) ]
	end
	if(!(constant_pert == 0.0))
		temp .+= constant_pert
	end


	for (iff,field) in enumerate(fields)
		m=getfield(mod, field)
		for i in eachindex(m)
			if(((m[i]-onlyin[iff][1])*(m[i]-onlyin[iff][2]))<0.0)
				m[i] += temp[i]
			end
		end
	end

	# random noise
	for field in fields
		f0 = Symbol((replace("$(field)", "χ" => "")),0)
		m=getfield(mod, field)
		m0=getfield(mod, f0)
		for i in eachindex(m)
			m[i] += randn() * randn_perc * 1e-2
		end

	end
	return mod
end


"""
Apply smoothing to `Seismic` using a Gaussian filter of zwidth and xwidth

# Arguments

* `mod::Seismic` : argument that is modified
* `zperc::Real` : smoothing percentage in z-direction
* `xperc::Real=zperc` : smoothing percentage in x-direction

# Keyword Arguments

* `zmin::Real=mod.mgrid[1][1]` : 
* `zmax::Real=mod.mgrid[1][end]` : 
* `xmin::Real=mod.mgrid[2][1]` : 
* `xmax::Real=mod.mgrid[2][end]` : 
* `fields` : fields of seismic model that are to be smooth
"""
function Seismic_smooth(mod::Seismic, zperc::Real, xperc::Real=zperc;
		 zmin::Real=mod.mgrid[1][1], zmax::Real=mod.mgrid[1][end],
		 xmin::Real=mod.mgrid[2][1], xmax::Real=mod.mgrid[2][end],
		 fields=[:χvp, :χρ, :χvs]
			)
	warn("ImageFiltering bug, not active")
	xwidth = Float64(xperc) * 0.01 * abs(mod.mgrid[2][end]-mod.mgrid[2][1])
	zwidth = Float64(zperc) * 0.01 * abs(mod.mgrid[1][end]-mod.mgrid[1][1])
	xnwin=Int(div(xwidth,step(mod.mgrid[2])*2.))
	znwin=Int(div(zwidth,step(mod.mgrid[1])*2.))

	izmin = Interpolation.indminn(mod.mgrid[1], Float64(zmin), 1)[1]
	izmax = Interpolation.indminn(mod.mgrid[1], Float64(zmax), 1)[1]
	ixmin = Interpolation.indminn(mod.mgrid[2], Float64(xmin), 1)[1]
	ixmax = Interpolation.indminn(mod.mgrid[2], Float64(xmax), 1)[1]

	modg=deepcopy(mod)
	for (i,iff) in enumerate(fields)
		m=view(getfield(mod, iff),izmin:izmax,ixmin:ixmax)
		mg=view(getfield(modg, iff),izmin:izmax,ixmin:ixmax)
#		imfilter!(mg, m, Kernel.gaussian([znwin,xnwin]));
	end
	return modg
end

"""
Return a truncated seismic model using input bounds.
Note that there is no interpolation going on here, but only truncation, so 
the input bounds cannot be strictly imposed.

# Arguments
* `mod::Seismic` : model that is truncated

# Keyword Arguments

* `zmin::Float64=mod.mgrid[1][1]` : 
* `zmax::Float64=mod.mgrid[1][end]` : 
* `xmin::Float64=mod.mgrid[2][1]` : 
* `xmax::Float64=mod.mgrid[2][end]` : 
"""
function Seismic_trun(mod::Seismic;
			 zmin::Float64=mod.mgrid[1][1], zmax::Float64=mod.mgrid[1][end],
			 xmin::Float64=mod.mgrid[2][1], xmax::Float64=mod.mgrid[2][end],
			 )

	izmin = Interpolation.indminn(mod.mgrid[1], zmin, 1)[1]
	izmax = Interpolation.indminn(mod.mgrid[1], zmax, 1)[1]
	ixmin = Interpolation.indminn(mod.mgrid[2], xmin, 1)[1]
	ixmax = Interpolation.indminn(mod.mgrid[2], xmax, 1)[1]

	x = mod.mgrid[2][ixmin:ixmax]
	z = mod.mgrid[1][izmin:izmax]

	# allocate model
	mod_trun = Seismic_zeros([z,x])
	adjust_bounds!(mod_trun, mod)

	for f in [:χvp, :χvs, :χρ]
		f0 = Symbol((replace("$(f)", "χ" => "")),0)
		setfield!(mod_trun, f, getfield(mod, f)[izmin:izmax, ixmin:ixmax])
		setfield!(mod_trun, f0, getfield(mod, f0)) # adjust bounds later after truncation if necessary
	end
	return mod_trun

end

"""
Extend a seismic model into PML layers
"""
function Seismic_pml_pad_trun(mod::Seismic, nlayer_rand, npml)

	mg=mod.mgrid
	mgex=[
        range(mg[1][1] - npml*step(mg[1]),
	     stop=mg[1][end] + npml*step(mg[1]), length=length(mg[1])+2*npml),
        range(mg[2][1] - npml*step(mg[2]),
	     stop=mg[2][end] + npml*step(mg[2]), length=length(mg[2])+2*npml)
	]

	modex=Seismic_zeros(mgex)
	adjust_bounds!(modex,mod)
	Seismic_pml_pad_trun!(modex, mod, nlayer_rand)

	return modex
end

"only padding implemented"
function Seismic_pml_pad_trun!(modex::Seismic, mod::Seismic, nlayer_rand)
	vp=mod.χvp; χ!(vp,mod.ref.vp,-1)
	vs=mod.χvs; χ!(vs,mod.ref.vs,-1)
	ρ=mod.χρ; χ!(ρ,mod.ref.ρ,-1)

	vpex=modex.χvp; vsex=modex.χvs; ρex=modex.χρ;
	pml_pad_trun!(vpex,vp,mod.vp0,nlayer_rand,50.0);	
	pml_pad_trun!(vsex,vs,[1,2],nlayer_rand,0.0);	
	pml_pad_trun!(ρex,ρ,mod.ρ0,nlayer_rand,0.0)
	χ!(vpex,modex.ref.vp,1); χ!(vsex,modex.ref.vs,1); χ!(ρex,modex.ref.ρ,1)
	χ!(vp,mod.ref.vp,1); χ!(vs,mod.ref.vs,1); χ!(ρ,mod.ref.ρ,1)
end

function pml_pad_trun(mod::Array{Float64,2}, np::Int64, flag::Int64=1)
	if(isequal(flag,1)) 
		nz = size(mod,1); nx = size(mod,2)
		modo = zeros(nz + 2*np, nx + 2*np)
		pml_pad_trun!(modo,mod,1)
	elseif(isequal(flag,-1))
		nz = size(mod,1); nx = size(mod,2)
		modo = zeros(nz - 2*np, nx - 2*np)
		pml_pad_trun!(mod,modo,-1)
	end
	return modo
end


"""
PML Extend a model on all four sides
"""
function pml_pad_trun!(modex::Array{Float64,2}, mod::Array{Float64,2}, bounds, nlayer_rand, var_fact=1.)
	np=(size(modex,1)-size(mod,1) == size(modex,2)-size(mod,2)) ? 
			div(size(modex,2)-size(mod,2),2) : error("modex size")
	nz = size(mod,1); nx = size(mod,2)
	#=
	# first fill with random numbers
	Random.seed!(1)
	#for i in eachindex(modex)
	#	modex[i]=rand(Distributions.Uniform(bounds...))
	#end
	# experiments with random boundary conditions
	sample = v->rand(Distributions.Truncated(Distributions.Normal(mean(bounds),v),bounds[1],bounds[2]))
	for ix in 1+np:np+nx
		for iz in 1:np
			v=var_fact*(np-iz)/np+1e-5
			modex[iz,ix] = sample(v)
		end
		for iz in nz+1+np:size(modex,1)
			v=var_fact*(iz-np-nz-1)/np+1e-5
			modex[iz,ix] = sample(v)
		end
	end
	for iz in 1+np:np+nz
		for ix in 1:np
			v=var_fact*(np-ix)/np+1e-5
			modex[iz,ix] = sample(v)
		end
		for ix in nx+1+np:size(modex,2)
			v=var_fact*(ix-np-nz-1)/np+1e-5
			modex[iz,ix] = sample(v)
		end
	end

	for ix in 1:np
		for iz in 1:np
			v=var_fact*sqrt((np-iz)^2+(np-ix)^2)/np+1e-5
			modex[iz,ix] = sample(v)
		end
		for iz in nz+1+np:size(modex,1)
			v=var_fact*sqrt((iz-np-nz-1)^2+(np-ix)^2)/np+1e-5
			modex[iz,ix] = sample(v)
		end
	end
	for ix in np+nx+1:size(modex,1)
		for iz in 1:np
			v=var_fact*sqrt((np-iz)^2+(ix-np-nx-1)^2)/np+1e-5
			modex[iz,ix] = sample(v)
		end
		for iz in nz+1+np:size(modex,1)
			v=var_fact*sqrt((iz-np-nz-1)^2+(ix-np-nx-1)^2)/np+1e-5
			modex[iz,ix] = sample(v)
		end
	end

	=#
	for ix in 1:nx
		for iz in 1:nz
			modex[np+iz,np+ix]=mod[iz,ix]
		end
	end

	# flood borders
	for ix in 1:size(modex,2)
		for iz in 1:np-nlayer_rand
			modex[iz,ix] = modex[np+1-nlayer_rand,ix]
		end
		for iz in nz+1+np+nlayer_rand:size(modex,1)
			modex[iz,ix] = modex[nz+np+nlayer_rand,ix];
		end
	end
	for iz in 1:size(modex,1)
		for ix in 1:np-nlayer_rand
			modex[iz,ix] = modex[iz,np+1-nlayer_rand]
		end
		for ix in nx+1+np+nlayer_rand:size(modex,2)
			modex[iz,ix] = modex[iz,nx+np+nlayer_rand]
		end
	end
	return modex
end


"""
function to resample in the model domain

# Arguments
* `mod::Seismic` : model
* `modi::Seismic` : model after interpolation
"""
function interp_spray!(mod::Seismic, modi::Seismic, attrib::Symbol, Battrib::Symbol=:B2; pa=nothing)
	if(pa===nothing)
		pa=Interpolation.Kernel([mod.mgrid[2], mod.mgrid[1]], [modi.mgrid[2], modi.mgrid[1]], Battrib)
	end

	"loop over fields in `Seismic`"
	Interpolation.interp_spray!(mod.χvp, modi.χvp, pa, attrib)
	Interpolation.interp_spray!(mod.χvs, modi.χvs, pa, attrib)
	Interpolation.interp_spray!(mod.χρ, modi.χρ,  pa, attrib)

end

#
#macro inter

function save(mod::Seismic, folder; N=100)
	!(isdir(folder)) && error("invalid directory")
	error("need to be updated")

	nx=length(mod.mgrid[2])
	nz=length(mod.mgrid[1])
	n=max(nx,nz);
	fact=(n>N) ? round(Int,n/N) : 1
	#mgrid=resamp(mod.mgrid, step(mod.mgrid[2])*fact, step(mod.mgrid[1])*fact)
	x=mgrid[2]
	z=mgrid[1]
	nx=length(x)
	nz=length(z)
	modo=Seismic_zeros(mgrid)
	adjust_bounds!(modo, mod)
	interp_spray!(mod, modo, :interp)

	for m in [:vp, :ρ, :Zp]
		# save original gf
		file=joinpath(folder, string("im", m, ".csv"))
		CSV.write(file,DataFrame(hcat(repeat(z,outer=nx),
					      repeat(x,inner=nz),vec(Seismic_get(modo, m)))))
	end
end
	

end # module
