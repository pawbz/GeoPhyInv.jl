#module Media

"""
Mutable type for storing medium parameters.
```julia
mod=Medium(mgrid, names=[:vp,:rho])
```
This initializes a subsurface model with `:vp` and `:rho` parameters. Print the
names of the medium parameters stored in `mod`.
```julia
names(mod)
```

## Indexing
* `mod.mgrid` : returns the spatial-grid bundle
* `mod[:vp]` : P-wave velocity
* `mod[:rho]` : mass density
* `mod[:Zp]` : P-wave impedance 
* `mod[:K]` : bulk modulus (unrelaxed when considering attenuation) 
* `mod[:Q]` : quality factor (relaxation times are optimized to be constant over all frequencies; see Robertsson, et. al, 1994)
* `mod.ref` : reference medium parameters 
* `mod.bounds` : bounds of medium parameters

"""
mutable struct Medium{N} # N is 2 for 2D, and 3 for 3D
	mgrid::Vector{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}}
	m::NamedStack{Array{Float64,N}}
	bounds::NamedArrays.NamedArray{Array{Float64,1},1,Array{Array{Float64,1},1},Tuple{OrderedCollections.OrderedDict{Symbol,Int64}}}
	ref::NamedArray{Float64,1,Array{Float64,1},Tuple{OrderedCollections.OrderedDict{Symbol,Int64}}}
	# store some floats specific to the medium
	fc::NamedArrays.NamedArray{Float64,1,Array{Float64,1},Tuple{OrderedCollections.OrderedDict{Symbol,Int64}}}
	# store some integer constants
	ic::NamedArrays.NamedArray{Int64,1,Array{Int64,1},Tuple{OrderedCollections.OrderedDict{Symbol,Int64}}}
	# store 3D arrays, if any
	m3::NamedStack{Array{Float64,3}}
end

# default construction without fcs and ics, and m3
function Medium(mgrid,m::NamedStack{Array{Float64,N}},bounds,ref) where N
	return Medium(mgrid,m,bounds,ref,NamedArray([0.0],([:o],)),NamedArray([0],([:o],)),NamedArray([zeros(1,1,1)], ([:o],))) 
end

function NamedArrays.names(mod::Medium)
	return names(mod.m)
end
	

include("attenuation.jl")
include("base.jl")
include("gallery.jl")



#=


"""
Data type to represent a seismic model.
A contrast function for a model m is given by ``χ(m) = \frac{m}{m0}-1``.

# Fields

* `vp0::Vector{Float64}` : [vpmin, vpmax]
* `vs0::Vector{Float64}` : [vsmin, vsmax]
* `rho0::Vector{Float64}` : [rhomin, rhomax]
* `χvp::Array{Float64,2}` : two-dimensional contrast model (χ) for vp, for e.g., zeros(length(mgrid[1]), length(mgrid[2]))
* `χvs::Array{Float64}` : two-dimensional contrast model (χ) for vs, for e.g., zeros(length(mgrid[1]), length(mgrid[2]))
* `χrho::Array{Float64}` : two-dimensional contrast model (χ) for density, for e.g., zeros(length(mgrid[1]), length(mgrid[2]))
* `mgrid` : array of ranges to determine the dimensions of models
"""
#
#	vp0::Vector{Float64}
#	vs0::Vector{Float64}
#	rho0::Vector{Float64}
#	χvp::Array{Float64,2}
#	χvs::Array{Float64,2}
#	χrho::Array{Float64,2}
#	# all derived fields from previous
#	K0::Vector{Float64}
#	mu0::Vector{Float64}
#	KI0::Vector{Float64}
#	muI0::Vector{Float64}
#	rhoI0::Vector{Float64}
#	"adding conditions that are to be false while construction"
#	Medium(vp0,vs0,rho0,χvp,χvs,χrho,mgrid,K0,mu0,KI0,muI0,rhoI0,ref) = 
#		any([
#		     any(vp0.<0.0), 
#		     any(vs0.<0.0),
#		     any(rho0.<0.0),
#		     #all([all(vp0 .≠ 0.0), any(χ(χvp,vp0,-1) .< vp0[1])]), # check vp bounds
#		     #all([all(vp0 .≠ 0.0), any(χ(χvp,vp0,-1) .> vp0[2])]), # check vp bounds
#		     #all([all(vs0 .≠ 0.0), any(χ(χvs,vs0,-1) .< vs0[1])]), # check vs bounds
#		     #all([all(vs0 .≠ 0.0), any(χ(χvs,vs0,-1) .> vs0[2])]), # check vs bounds
#		     #all([all(rho0 .≠ 0.0), any(χ(χrho,rho0,-1) .< rho0[1])]), # check rho bounds
#		     #all([all(rho0 .≠ 0.0), any(χ(χrho,rho0,-1) .> rho0[2])]), # check rho bounds
#		     size(χvp) != (length(mgrid[1]), length(mgrid[2])), # dimension check
#		     size(χvs) != (length(mgrid[1]), length(mgrid[2])), # dimension check
#		     size(χrho) != (length(mgrid[1]), length(mgrid[2])) # dimension check
#		    ]) ? 
#		       error("error in Medium construction") : new(vp0,vs0,rho0,χvp,χvs,χrho,mgrid,K0,mu0,KI0,muI0,rhoI0,ref)
#end
#
# check whether a seismic model is bounded
function isbounded(mod::Medium)
	@warn "need to be implemented"
	result=true
end

# method to create Medium object without using derived fields
function Medium(vp0, vs0, rho0, χvp, χvs, χrho, 
		 mgrid::Vector{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}})
	mod=Medium(vp0, vs0, rho0, χvp, χvs, χrho, mgrid,
	    zeros(2), zeros(2), zeros(2), zeros(2), zeros(2),
	    Medium_ref(zeros(8)...))
	update_derived!(mod)
	return mod
end

=#

# update the derived fields of the Seisimc model
function update!(mod::Medium)
	"note that: mean(K0) ≠ 1/mean(KI0)"
	mb=deepcopy(mod.bounds)

	K = mb[:vp] .* mb[:vp] .* mb[:rho]; 
	KI = reverse(inv.(K),dims=1);
	rho = mb[:rho]; 
	rhoI = reverse(inv.(rho),dims=1);
	newnames=[:K,:KI,:rhoI]

	if(:vs ∈ names(mb,1))
		mu = mb[:vs] .* mb[:vs] .* mb[:rho];
		muI = reverse(inv.(mu),dims=1);
		mbnew=NamedArray([K,KI,rhoI,mu,muI],(vcat(newnames,[:mu,:muI]),))
	else
		mbnew=NamedArray([K,KI,rhoI],(newnames,))
	end


	for s in newnames
		# println(typeof(mod.bounds), '\'NamedArray([[0.0,0.0]],[[s],]))
		(s ∉ names(mod.bounds)[1]) && (mod.bounds=vcat(mod.bounds,NamedArray([[0.0,0.0]],[s,])))
	end

	for s in newnames
		copyto!(mod.bounds[s],mbnew[s])
	end

	for s in names(mod.bounds)[1]
		meanvalue=Statistics.mean(mod.bounds[s])
		if(s ∈ names(mod.ref)[1])
			mod.ref[s]=meanvalue
		else
			mod.ref=vcat(mod.ref, NamedArray([meanvalue], [s,]))
		end
	end

	return nothing
end


"""
Return medium property bounds based on maximum and minimum values of the array and frac.
The bounds cannot be less than zero
"""
function bounds(mod::AbstractArray, frac::Real=0.1)
	any(mod .< 0.0) && error("model values less than zero")
	bounds=zeros(2)
	bound=frac*Statistics.mean(mod)
	bounds[1] = ((minimum(mod) - bound)<0.0) ? 0.0 : (minimum(mod) - bound)
	bounds[2] = maximum(mod)+bound
	return bounds
end


update!(mod::Medium, mod0::Medium)=update!(mod, names(mod0.bounds)[1], mod0.bounds)

function update!(mod::Medium, fields::AbstractVector{Symbol}, bounds)
	@assert !(bounds===nothing)
	for name in names(mod.bounds)[1]
		if(name ∈ fields)
			i=findall(x->x==name, fields)[1]
			copyto!(mod.bounds[name], Float64.(bounds[i]))
		end
	end
	update!(mod)
end

"""
Adjust the bounds and hence the reference values.
Since the reference values are adjust the χ fields should also be changed
"""
function update!(mod::Medium, frac::Real)# 0.1
	for name in names(mod.m)[1]
		b = bounds(mod.m[name], frac)
		copyto!(mod.bounds[name], Float64.(b))
	end
	update!(mod)
	return mod
end


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
In-place method to add features to the input model.

# Arguments
* `mod::Medium` : model that is modified

# Keyword Arguments

* `point_loc::Vector{Float64}=[0., 0.,]` : approx location of point pert
* `point_pert::Float64=0.0` : perturbation at the point scatterer
* `ellip_loc::Vector{Float64}=nothing` : location of center of perturbation, [z, x]
* `ellip_rad::Float64=0.0` : size of elliptic perturbation
* `ellip_pert::Float64=0.1` : perturbation inside the ellipse
* `ellip_α=0.0` : rotate the ellipse
* `rect_loc::Array{Float64}=nothing` : rectangle location, [zmin, xmin, zmax, xmax]
* `rect_pert::Float64=0.1` : perturbation in a rectangle
* `constant_pert::Float64=0.0` : constant perturbation 
* `randn_pert::Float64=0.0` : percentage of reference values for additive random noise
* `fields::Vector{Symbol}=[:χvp,:χrho,:χvs]` : which fields are to be modified?
* `onlyin` : `mod` is modified only when field values are in these ranges 
"""
function update!(mod::Medium, fields::Vector{Symbol}; 
		       point_loc=fill(0.0,ndims(mod)),
		       point_pert::Float64=0.0,
		       ellip_loc=fill(0.0,ndims(mod)),
		       ellip_rad=0.0,
		       ellip_pert::Real=0.0,
		       ellip_α=0.0,
		       rect_loc=fill(0.0,ndims(mod)),
		       rect_pert::Float64=0.0,
		       constant_pert::Float64=0.0,
		       randn_perc::Real=0.0,
			   layerlocations=nothing
			   )

    # only editing basic fields is allowed
	for field in fields
		@assert field ∈ [:vp,:vs,:rho,:Q]
	end
	rect_loc=convert.(Float64,rect_loc);
	ellip_loc=convert.(Float64,ellip_loc);

	temp = zeros(length.(mod.mgrid)...)

	ipointloc = [Interpolation.indminn(mod.mgrid[i], Float64(point_loc[i]), 1)[1] for i in 1:length(mod.mgrid)]
	temp[ipointloc...] += point_pert

	if(!(ellip_pert == 0.0))
		α=ellip_α*pi/180.
		# circle or ellipse
		rads= (length(ellip_rad)==1) ? [ellip_rad[1],ellip_rad[1]] : [ellip_rad[1],ellip_rad[2]]

		temp += [(((((mod.mgrid[2][ix]-ellip_loc[2])*cos(α)+(mod.mgrid[1][iz]-ellip_loc[1])*sin(α))^2*inv(rads[1]^2) + 
	      ((-mod.mgrid[1][iz]+ellip_loc[1])*cos(α)+(mod.mgrid[2][ix]-ellip_loc[2])*sin(α))^2*inv(rads[2]^2)) <= 1.) ? Float64(ellip_pert) : 0.0)  for 
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
		m=mod.m[field]
		if(!(layerlocations===nothing))
			onlyatvalues=[]
			for ipos in layerlocations
				ipx = Interpolation.indminn(mod.mgrid[2], Float64(ipos[2]), 1)[1] 
				ipz = Interpolation.indminn(mod.mgrid[1], Float64(ipos[1]), 1)[1] 
				push!(onlyatvalues,m[ipz,ipx])
			end
		end
		for i in eachindex(m)
			# if(((m[i]-onlyin[iff][1])*(m[i]-onlyin[iff][2]))<0.0)
			if(layerlocations===nothing)
				m[i] += temp[i]
			else
				if( m[i] ∈ onlyatvalues)
					m[i] += temp[i]
				end
			end
		end
	end

	# random noise
	for field in fields
		m=mod.m[field]
		m0=mod.ref[field]
		for i in eachindex(m)
			m[i] = χ(χ(m[i],m0,1) + randn() * randn_perc * 1e-2, m0, -1)
		end

	end
	return nothing
end

include("chainrule.jl")



"""
Apply smoothing to `Medium` using a Gaussian filter of zwidth and xwidth

# Arguments

* `mod::Medium` : argument that is modified
* `zperc::Real` : smoothing percentage in z-direction
* `xperc::Real=zperc` : smoothing percentage in x-direction

# Keyword Arguments

* `zmin::Real=mod.mgrid[1][1]` : 
* `zmax::Real=mod.mgrid[1][end]` : 
* `xmin::Real=mod.mgrid[2][1]` : 
* `xmax::Real=mod.mgrid[2][end]` : 
* `fields` : fields of seismic model that are to be smooth
"""
function Medium_smooth(mod::Medium, zperc::Real, xperc::Real=zperc;
		 zmin::Real=mod.mgrid[1][1], zmax::Real=mod.mgrid[1][end],
		 xmin::Real=mod.mgrid[2][1], xmax::Real=mod.mgrid[2][end],
		 fields=[:vp, :rho]
			)
	xwidth = Float64(xperc) * 0.01 * abs(mod.mgrid[2][end]-mod.mgrid[2][1])
	zwidth = Float64(zperc) * 0.01 * abs(mod.mgrid[1][end]-mod.mgrid[1][1])
	xnwin=Int(div(xwidth,step(mod.mgrid[2])*2.))
	znwin=Int(div(zwidth,step(mod.mgrid[1])*2.))

	izmin = Interpolation.indminn(mod.mgrid[1], Float64(zmin), 1)[1]
	izmax = Interpolation.indminn(mod.mgrid[1], Float64(zmax), 1)[1]
	ixmin = Interpolation.indminn(mod.mgrid[2], Float64(xmin), 1)[1]
	ixmax = Interpolation.indminn(mod.mgrid[2], Float64(xmax), 1)[1]
	
	# @warn "check this routine, smooth contrast values instead?"

	modg=deepcopy(mod)
	for (i,iff) in enumerate(fields)
		m=view(mod[iff],izmin:izmax,ixmin:ixmax)
		mg=view(modg[iff],izmin:izmax,ixmin:ixmax)
		imfilter!(mg, m, Kernel.gaussian([znwin,xnwin]));
	end
	return modg
end

"""
Return a truncated seismic model using input bounds.
Note that there is no interpolation going on here, but only truncation, so 
the input bounds cannot be strictly imposed.

# Arguments
* `mod::Medium` : model that is truncated

# Keyword Arguments

* `zmin::Real=mod.mgrid[1][1]` : 
* `zmax::Real=mod.mgrid[1][end]` : 
* `xmin::Real=mod.mgrid[2][1]` : 
* `xmax::Real=mod.mgrid[2][end]` : 
"""
function Medium_trun(mod::Medium;
			 zmin::Real=mod.mgrid[1][1], zmax::Real=mod.mgrid[1][end],
			 xmin::Real=mod.mgrid[2][1], xmax::Real=mod.mgrid[2][end],
			 )

	izmin = Interpolation.indminn(mod.mgrid[1], Float64(zmin), 1)[1]
	izmax = Interpolation.indminn(mod.mgrid[1], Float64(zmax), 1)[1]
	ixmin = Interpolation.indminn(mod.mgrid[2], Float64(xmin), 1)[1]
	ixmax = Interpolation.indminn(mod.mgrid[2], Float64(xmax), 1)[1]

	x = mod.mgrid[2][ixmin:ixmax]
	z = mod.mgrid[1][izmin:izmax]

	# allocate model
	mod_trun = Medium([z,x], names(mod.m)[1])
	update!(mod_trun, mod)

	for f in names(mod.m)[1]
		m=view(mod[f],izmin:izmax, ixmin:ixmax)
		copyto!(mod_trun[f],m)
	end
	return mod_trun

end


"""
Extend a seismic model into PML layers
"""
function Medium_pml_pad_trun(mod::Medium, nlayer_rand, npml)

	mg=mod.mgrid
	mgex=[
        range(mg[1][1] - npml*step(mg[1]),
	     stop=mg[1][end] + npml*step(mg[1]), length=length(mg[1])+2*npml),
        range(mg[2][1] - npml*step(mg[2]),
	     stop=mg[2][end] + npml*step(mg[2]), length=length(mg[2])+2*npml)
	]

	modex=Medium(mgex, names(mod.m)[1])
	update!(modex,mod)
	Medium_pml_pad_trun!(modex, mod, nlayer_rand)

	return modex
end

"only padding implemented"
function Medium_pml_pad_trun!(modex::Medium, mod::Medium, nlayer_rand)
	vp=mod[:vp];
	rho=mod[:rho]; 

	vpex=modex[:vp]; rhoex=modex[:rho];
	pml_pad_trun!(vpex,vp,mod.bounds[:vp],nlayer_rand,50.0);	
	pml_pad_trun!(rhoex,rho,mod.bounds[:rho],nlayer_rand,0.0)
	if(:vs ∈ names(mod.m)[1])
		vs=mod[:vs];
		vsex=modex[:vs]; 
		pml_pad_trun!(vsex,vs,mod.bounds[:vs],nlayer_rand,50.0);	
	end
	if(:Q ∈ names(mod.m)[1])
		Q=mod[:Q];
		Qex=modex[:Q]; 
		pml_pad_trun!(Qex,Q,mod.bounds[:Q],nlayer_rand,50.0);	
	end
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
* `mod::Medium` : model
* `modi::Medium` : model after interpolation
"""
function interp_spray!(mod::Medium, modi::Medium, attrib::Symbol, Battrib::Symbol=:B2, fields=[:vp,:rho]; pa=nothing)
	if(pa===nothing)
		pa=Interpolation.Kernel([mod.mgrid[2], mod.mgrid[1]], [modi.mgrid[2], modi.mgrid[1]], Battrib)
	end

	# :vs is missing here....
	"loop over fields in `Medium`, add vs later"
	for field in fields
		Interpolation.interp_spray!(mod[field],  modi[field], pa, attrib)
	end
end

#
#macro inter

function save(mod::Medium, folder; N=100)
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
	modo=Medium_zeros(mgrid)
	adjust_bounds!(modo, mod)
	interp_spray!(mod, modo, :interp)

	for m in [:vp, :rho, :Zp]
		# save original gf
		file=joinpath(folder, string("im", m, ".csv"))
		CSV.write(file,DataFrame(hcat(repeat(z,outer=nx),
				repeat(x,inner=nz),vec(modo[m]))))
	end
end
	
