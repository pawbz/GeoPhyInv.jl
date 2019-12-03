"""
Initialize an model with zeros, given mgrid and names.
"""
function Base.zero(::Type{Medium}, mgrid, names=[:vp,:rho])
	nz=length(mgrid[1]); nx=length(mgrid[2])
	return Medium(mgrid,NamedArray([zeros(nz,nx) for i in names], (names,)),
	               zeros(nz,nx), 
		       NamedArray([fill(0.0,2) for i in names], (names,)), 
		       NamedArray([0.0 for i in names],(names,)))
end


"""
Return true if a `Seismic` model is just allocated with zeros.
"""
function Base.iszero(mod::Medium)
	return any((mod.vp0 .* mod.rho0) .== 0.0) ? true : false
end


"""
Return a similar model to the input model, used for allocation.
"""
function Base.similar(mod::Medium)
	return zero(Medium, mod.mgrid, names(mod.m)[1])
end


"Compare if two `Seismic` models are equal"
function Base.isequal(mod1::T, mod2::T) where {T<:Union{Medium}}

	fnames = fieldnames(T)
	vec=[(isequal(getfield(mod1, name),getfield(mod2, name))) for name in fnames]
	return all(vec)
end

"""
Return if two `Seismic` models have same dimensions and bounds.
"""
function Base.isapprox(mod1::Medium, mod2::Medium)
	return isequal(mod1.mgrid,mod2.mgrid) & isequal(mod1.ref,mod2.ref) & isequal(size.(mod1.m), size.(mod2.m))
end

"""
Copy for `Seismic` models. The models should have same bounds and sizes.
Doesn't allocate any memory.
"""
function Base.copyto!(modo::Medium, mod::Medium)
	@assert isapprox(modo, mod)
	for name in names(mod.m)[1]
		copyto!(modo.m[name], mod.m[name])
	end
	return modo
end


#function Base.fill!(mod::Medium, k::Float64=0.0)
#	for m in mod.χ
#		fill!(m,k)
#	end
#	return mod
#end


"""
Print information about `Seismic`
"""
function Base.print(mod::Medium, name::String="")
	println("\tSeismic Medium:\t",name)
	println("\t> number of samples:\t","x\t",length(mod.mgrid[2]),"\tz\t",length(mod.mgrid[1]))
	println("\t> sampling intervals:\t","x\t",step(mod.mgrid[2]),"\tz\t",step(mod.mgrid[1]))
	println("\t> vp:\t","min\t",minimum(Seismic_get(mod,:vp)),"\tmax\t",maximum(Seismic_get(mod,:vp)))
	println("\t> vp bounds:\t","min\t",mod.vp0[1],"\tmax\t",mod.vp0[2])
	println("\t> rho:\t","min\t",minimum(Seismic_get(mod,:rho)),"\tmax\t",maximum(Seismic_get(mod,:rho)))
	println("\t> rho bounds:\t","min\t",mod.rho0[1],"\tmax\t",mod.rho0[2])
end

"""
Return `Seismic` with zeros everywhere;
this method is used for preallocation.

# Arguments
* `mgrid` : used for sizes of χ fields 
"""
function Seismic_zeros(mgrid::Vector{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}})
	nz=length(mgrid[1]); nx=length(mgrid[2])
	return Seismic(mgrid,NamedArray([zeros(nz,nx),zeros(nx,nz)], ([:vp,:rho],)),
		       NamedArray([fill(0.0,2),fill(0.0,2)], ([:vp,:rho],)), 
		       NamedArray([0.0,0.0],([:vp,:rho],)))
end

