__precompile__()

module Models


import JuMIT.Grid
import JuMIT.Interpolation
import JuMIT.Smooth
using DataFrames
using CSV

"""
Data type fo represent a seismic model.
A contrast function for a model m is given by ``χ(m) = \frac{m}{m0}-1``.

# Fields

* `vp0::Vector{Float64}` : [vpmin, vpmax]
* `vs0::Vector{Float64}` : [vsmin, vsmax]
* `ρ0::Vector{Float64}` : [ρmin, ρmax]
* `χvp::Array{Float64,2}` : two-dimensional contrast model (χ) for vp, for e.g., zeros(mgrid.nz, mgrid.nx)
* `χvs::Array{Float64}` : two-dimensional contrast model (χ) for vs, for e.g., zeros(mgrid.nz, mgrid.nx)
* `χρ::Array{Float64}` : two-dimensional contrast model (χ) for density, for e.g., zeros(mgrid.nz, mgrid.nx)
* `mgrid::Grid.M2D` : two-dimensional grid to determine the dimensions of models
"""
type Seismic
	vp0::Vector{Float64}
	vs0::Vector{Float64}
	ρ0::Vector{Float64}
	χvp::Array{Float64,2}
	χvs::Array{Float64,2}
	χρ::Array{Float64,2}
	mgrid::Grid.M2D
	"adding conditions that are to be false while construction"
	Seismic(vp0,vs0,ρ0,χvp,χvs,χρ,mgrid) = 
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
		     size(χvp) != (length(mgrid.z), length(mgrid.x)), # dimension check
		     size(χvs) != (length(mgrid.z), length(mgrid.x)), # dimension check
		     size(χρ) != (length(mgrid.z), length(mgrid.x)) # dimension check
		    ]) ? 
		       error("error in Seismic construction") : new(vp0,vs0,ρ0,χvp,χvs,χρ,mgrid)
end

"""
Print information about `Seismic`
"""
function Base.print(mod::Seismic, name::String="")
	println("\tSeismic Model:\t",name)
	println("\t> number of samples:\t","x\t",mod.mgrid.nx,"\tz\t",mod.mgrid.nz)
	println("\t> sampling intervals:\t","x\t",mod.mgrid.δx,"\tz\t",mod.mgrid.δz)
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
	bound=frac*mean(mod)
	bounds[1] = ((minimum(mod) - bound)<0.0) ? 0.0 : (minimum(mod) - bound)
	bounds[2] = maximum(mod)+bound
	return bounds
end

"""
Adjust the bounds and hence the reference values.
Since the reference values are adjust the χ fields should also be changed
"""
function adjust_bounds!(mod::Seismic, frac::Float64=0.1)
	for f in [:χvp, :χρ, :χvs]
		f0 = Symbol((replace("$(f)", "χ", "")),0)
		fm = Symbol((replace("$(f)", "χ", "")))
		m = Seismic_get(mod, fm)
		m0 = bounds(m, frac)
		setfield!(mod, f0, m0)
		setfield!(mod, f, χ(m, m0, 1))
	end
	return mod
end

"""
Return `Seismic` with zeros everywhere;
this method is used for preallocation.

# Arguments
* `mgrid::Grid.M2D` : used for sizes of χ fields 
"""
function Seismic_zeros(mgrid::Grid.M2D)
	return Seismic(fill(0.0,2), fill(0.0,2), fill(0.0,2), 
		zeros(mgrid.nz, mgrid.nx), zeros(mgrid.nz, mgrid.nx),
		zeros(mgrid.nz, mgrid.nx), deepcopy(mgrid))
end

function Base.fill!(mod::Seismic, k::Float64=0.0)
	mod.χvp[:]=k
	mod.χρ[:]=k
	mod.χvs[:]=k
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
function Base.isequal(mod1::Seismic, mod2::Seismic)
	fnames = fieldnames(Seismic)
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
"""
function Base.copy!(modo::Seismic, mod::Seismic)
	if(isapprox(modo, mod))
		for f in [:χvp, :χρ, :χvs]
			modof=getfield(modo,f)
			modf=getfield(mod,f)
			for i in eachindex(modf)
				modof[i]=modf[i]
			end
		end
		return modo
	else
		error("attempt to copy dissimilar models")
	end
end

function Seismic_get(mod::Seismic, attrib::Symbol)
	K0 = mod.vp0 .* mod.vp0 .* mod.ρ0; 
	KI0 = flipdim(K0.^(-1.0),1);
	"note that: mean(K0) ≠ 1/mean(KI0)"
	ρ0 = mod.ρ0; ρI0 = flipdim(ρ0.^(-1.0),1);
	μ0 = mod.vs0 .* mod.vs0 .* mod.ρ0;
	if(attrib == :ρ0)
		return mod.ρ0;
	elseif(attrib == :vp0)
		return mod.vp0;
	elseif(attrib == :vs0)
		return mod.vs0;
	elseif(attrib == :ρI0)
		return ρI0
	elseif(attrib == :K0)
		return K0
	elseif(attrib == :μ0)
		return μ0
	elseif(attrib == :KI0)
		return KI0
	else
		# allocate
		rout=zeros(mod.mgrid.nz, mod.mgrid.nx)
		Seismic_get!(rout, mod, attrib)
		return rout
	end
end

"""
Get other dependent model parameters of a seismic model
that are not present in `Seismic`.

* `:ρI` : inverse of density
* `:Zp` : P-wave impedance
"""
function Seismic_get!(modmat::Array{Float64,2}, mod::Seismic, attrib::Symbol)
	K0 = mod.vp0 .* mod.vp0 .* mod.ρ0; 
	KI0 = flipdim(K0.^(-1.0),1);
	"note that: mean(K0) ≠ 1/mean(KI0)"
	ρ0 = mod.ρ0; ρI0 = flipdim(ρ0.^(-1.0),1);
	μ0 = mod.vs0 .* mod.vs0 .* mod.ρ0;
	ρ=mod.χρ; χ!(ρ, mod.ρ0,-1) # undo it later
	vp=mod.χvp; χ!(vp, mod.vp0,-1) # undo it later
	vs=mod.χvs; χ!(vs, mod.vs0,-1) # undo it later
	if(attrib == :ρI)
		for i in eachindex(modmat)
			modmat[i]=inv(ρ[i])
		end
	elseif(attrib == :ρ)
		copy!(modmat, ρ)
	elseif(attrib == :vp)
		copy!(modmat, vp)
	elseif(attrib == :vs)
		copy!(modmat, vs)
	elseif(attrib == :χρI)
		for i in eachindex(modmat)
			modmat[i]=inv(ρ[i])
		end
		χ!(modmat, ρI0, 1)
	elseif(attrib == :χρ)
		copy!(modmat, ρ)
		χ!(modmat, mod.ρ0, 1)
	elseif(attrib == :χvp)
		copy!(modmat, vp)
		χ!(modmat, mod.vp0, 1)
	elseif(attrib == :χK)
		for i in eachindex(modmat)
			modmat[i]= vp[i]*vp[i]*ρ[i]
		end
		χ!(modmat, mod.K0, 1)
	elseif(attrib == :χμ)
		for i in eachindex(modmat)
			modmat[i]= vs[i]*vs[i]*ρ[i]
		end
		χ!(modmat, mod.μ0, 1)
	elseif(attrib == :χKI)
		for i in eachindex(modmat)
			modmat[i]=inv(vp[i]*vp[i]*ρ[i])
		end
		χ!(modmat, KI0, 1)
	elseif(attrib == :KI)
		for i in eachindex(modmat)
			modmat[i]=inv(vp[i]*vp[i]*ρ[i])
		end
	elseif(attrib == :K)
		for i in eachindex(modmat)
			modmat[i]= vp[i]*vp[i]*ρ[i]
		end
	elseif(attrib == :Zp)
		for i in eachindex(modmat)
			modmat[i]= vp[i]*ρ[i]
		end
	else
		error("invalid attrib")
	end
	χ!(ρ, mod.ρ0,1) 
	χ!(vp, mod.vp0,1)
	χ!(vs, mod.vp0,1)
	return modmat
end

"""
Re-parameterization routine 
that modifies the fields 
`χvp` and `χρ` of an input seismic model
using two input vectors.

# Arguments

* `mod::Seismic` : to be updated
* `x1::Array{Float64,2}` : contrast of inverse bulk modulus
* `x2::Array{Float64,2}` : contrast of inverse density
* `attribs:::Vector{Symbol}` : [:χKI, :χρI]
"""
function Seismic_reparameterize!(
	              mod::Seismic,
		      x1::Array{Float64,2},
		      x2::Array{Float64,2},
		      attribs::Vector{Symbol}=[:χKI, :χρI]
		      )
	iszero(mod) ? error("mod cannot be zero") : nothing
	size(x1) == (mod.mgrid.nz, mod.mgrid.nx) ? nothing : error("size x1")
	size(x2) == (mod.mgrid.nz, mod.mgrid.nx) ? nothing : error("size x2")
	if(attribs == [:χKI, :χρI]) 
		ρ = (χ(x2, Seismic_get(mod, :ρI0), -1)).^(-1);
		K = (χ(x1, Seismic_get(mod, :KI0), -1)).^(-1);
		vp = sqrt.((K ./ ρ));
		mod.χvp = copy(χ(vp, mod.vp0, 1));
		mod.χρ = copy(χ(ρ, mod.ρ0, 1));
	else
		error("invalid attribs")
	end
	return mod
end

"""
Use chain rule to output gradients with 
respect to χvp and χρ from  gradients 
with respect to KI and ρI.

# Arguments

* `gmod::Seismic` : gradient model
* `mod::Seismic` : model required for chain rule
* `g1` : gradient of an objective function with respect `attribs[1]`
* `g2` : gradient of an objective function with respect `attribs[2]`
* `attribs::Vector{Symbol}=[:χKI, :χρI]` :  
* `flag::Int64=1` :
  * `=1` updates `gmod` using `g1` and `g2`
  * `=-1` updates `g1` and `g2` using `gmod`
"""
function Seismic_chainrule!(
		      gmod::Seismic,
		      mod::Seismic,
		      g1::AbstractVector{Float64},
		      g2::AbstractVector{Float64},
		      attribs::Vector{Symbol}=[:χKI, :χρI],
		      flag::Int64=1
		      )

	if(flag == 1)
		if(attribs == [:χKI, :χρI]) 
			vp0 = mod.vp0;	vs0 = mod.vs0;	ρ0 = mod.ρ0
			ρ=view(mod.χρ,:); χ!(ρ, mod.ρ0,-1) # undo it later
			vp=view(mod.χvp,:); χ!(vp, mod.vp0,-1) # undo it later
			KI0=Seismic_get(mod,:KI0)
			ρI0=Seismic_get(mod,:ρI0)
			χg!(g1,KI0,-1)
			χg!(g2,ρI0,-1)

			gχρ=view(gmod.χρ,:)
			gχvp=view(gmod.χvp,:)
			# gradient w.r.t  χρ
			@. gχρ = -1.*inv(ρ)^2*g2-inv(vp*vp*ρ*ρ)*g1
			# gradient w.r.t χvp
			@. gχvp = -2.*inv(vp^3)*inv(ρ)*g1
			χg!(gχρ,ρ0,1)
			χg!(gχvp,vp0,1)
			# are these necessary?
			#gmod.χvs = copy(zeros(gmod.χvp))
			#gmod.mgrid = deepcopy(mod.mgrid);	
			#gmod.vp0 = copy(mod.vp0);	gmod.ρ0 = copy(mod.ρ0)
			χ!(ρ, mod.ρ0,1)
			χ!(vp, mod.vp0,1)
			χg!(g1,KI0,1)
			χg!(g2,ρI0,1)
		else
			error("invalid attribs")
		end
	elseif(flag == -1)
		χg!(gmod.χvp, mod.vp0, -1) # undo it later
		χg!(gmod.χρ, mod.ρ0, -1) # undo it later
		gvp = view(gmod.χvp, :)
		gρ = view(gmod.χρ, :)
		K0 = mean(Seismic_get(mod, :K0))
		KI0 = mean(Seismic_get(mod, :KI0))
		ρ0 = mean(mod.ρ0)
		ρI0 = mean(Seismic_get(mod, :ρI0))
		ρ0 = mean(mod.ρ0)
		ρ=view(mod.χρ,:); χ!(ρ, mod.ρ0,-1) # undo it later
		vp=view(mod.χvp,:); χ!(vp, mod.vp0,-1) # undo it later

		nznx = mod.mgrid.nz * mod.mgrid.nx;
		if(attribs == [:χK, :χρ])
			@. g1 = gvp * inv(ρ) * inv(vp) 
			scale!(gmod.χvp, 0.5*K0)
			@. g2 = gp
			@. g2 = gvp * vp * inv(ρ)
			scale!(g2,-0.5)
			@. g2 += gρ
			scale!(g2,ρ0)
		elseif(attribs == [:χKI, :χρI])
			chain_g1!(g1,gvp,ρ,vp,KI0)
			chain_g2!(g2,gρ,ρI0,ρ,gvp,vp)
		else
			error("invalid attribs")
		end
		χg!(gvp,mod.vp0,1)
		χg!(gρ,mod.ρ0,1)
		χ!(ρ, mod.ρ0,1)
		χ!(vp, mod.vp0,1)
	else
		error("invalid flag")
	end

end

function chain_g1!(g1,gvp,ρ,vp,KI0)
	@simd for i in eachindex(g1)
		g1[i] = -0.5 * KI0 * gvp[i] * ρ[i] *  vp[i]^3
	end
end
function chain_g2!(g2,gρ,ρI0,ρ,gvp,vp)
	@simd for i in eachindex(g2)
		g2[i] = (gρ[i] * ρI0 * -1. * ρ[i]*ρ[i]) + 
			(0.5 * ρI0 .* gvp[i] .*  vp[i] .* ρ[i])
	end
end




"""
Return dimensionless contrast model parameter using the
reference value.

# Arguments
* `mod::Array{Float64}` : subsurface parameter
* `mod0::Vector{Float64}` : reference value is mean of this vector
* `flag::Int64=1` : 
"""
function χ(mod::AbstractArray{Float64}, mod0::Vector{Float64}, flag::Int64=1)
	mod_out=copy(mod)
	χ!(mod_out, mod0, flag)
	return mod_out
end
function χ!(mod::AbstractArray{Float64}, mod0::Vector{Float64}, flag::Int64=1)
	m0 = mean(mod0)
	if(flag == 1)
		if(iszero(m0))
			return nothing
		else
			@simd for i in eachindex(mod)
				mod[i] = ((mod[i] - m0) * inv(m0))
			end
			return mod 
		end
	elseif(flag == -1)
		@simd for i in eachindex(mod)
			mod[i] = mod[i] * m0 + m0
		end
		return mod 
	end
end


"""
Gradients
Return contrast model parameter using the
reference value.
"""
function χg(mod::AbstractArray{Float64}, mod0::Vector{Float64}, flag::Int64=1)
	mod_out=copy(mod)
	χg!(mod_out, mod0, flag)
	return mod_out
end
function χg!(mod::AbstractArray{Float64}, mod0::Vector{Float64}, flag::Int64=1)
	if(flag == 1)
		m0 = mean(mod0)
		scale!(mod, m0)
	elseif(flag == -1)
		m0 = inv(mean(mod0))
		scale!(mod, m0)
	end
	return mod
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
		       point_loc::Vector{Float64}=[0., 0.,],
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

	temp = zeros(mod.mgrid.nz, mod.mgrid.nx)

	ipointlocx = Interpolation.indminn(mod.mgrid.x, point_loc[2], 1)[1] 
	ipointlocz = Interpolation.indminn(mod.mgrid.z, point_loc[1], 1)[1] 
	temp[ipointlocz, ipointlocx] += point_pert

	if(!(ellip_pert == 0.0))
		α=ellip_α*pi/180.
		# circle or ellipse
		rads= (length(ellip_rad)==1) ? [ellip_rad[1],ellip_rad[1]] : [ellip_rad[1],ellip_rad[2]]

		temp += [(((((mod.mgrid.x[ix]-ellip_loc[2])*cos(α)+(mod.mgrid.z[iz]-ellip_loc[1])*sin(α))^2*inv(rads[1]^2) + 
	      ((-mod.mgrid.z[iz]+ellip_loc[1])*cos(α)+(mod.mgrid.x[ix]-ellip_loc[2])*sin(α))^2*inv(rads[2]^2)) <= 1.) ? ellip_pert : 0.0)  for 
	   		iz in 1:mod.mgrid.nz, ix in 1:mod.mgrid.nx ]
	end
	if(!(rect_pert == 0.0))
		temp += [
			(((mod.mgrid.x[ix]-rect_loc[4]) * (mod.mgrid.x[ix]-rect_loc[2]) < 0.0) & 
			((mod.mgrid.z[iz]-rect_loc[3]) * (mod.mgrid.z[iz]-rect_loc[1]) < 0.0)) ?
			rect_pert : 0.0  for
			iz in 1:mod.mgrid.nz, ix in 1:mod.mgrid.nx ]
	end
	if(!(constant_pert == 0.0))
		temp += constant_pert
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
		f0 = Symbol((replace("$(field)", "χ", "")),0)
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
* `zperc::Float64` : smoothing percentage in z-direction
* `xperc::Float64=zperc` : smoothing percentage in x-direction

# Keyword Arguments

* `zmin::Float64=mod.mgrid.z[1]` : 
* `zmax::Float64=mod.mgrid.z[end]` : 
* `xmin::Float64=mod.mgrid.x[1]` : 
* `xmax::Float64=mod.mgrid.x[end]` : 
"""
function Seismic_smooth!(mod::Seismic, zperc::Float64, xperc::Float64=zwidth;
			 zmin::Float64=mod.mgrid.z[1], zmax::Float64=mod.mgrid.z[end],
			 xmin::Float64=mod.mgrid.x[1], xmax::Float64=mod.mgrid.x[end] 
			)
	xwidth = xperc * 0.01 * abs(mod.mgrid.x[end]-mod.mgrid.x[1])
	zwidth = zperc * 0.01 * abs(mod.mgrid.z[end]-mod.mgrid.z[1])
	xnwin=Int(div(xwidth,mod.mgrid.δx*2.))
	znwin=Int(div(zwidth,mod.mgrid.δz*2.))

	izmin = Interpolation.indminn(mod.mgrid.z, zmin, 1)[1]
	izmax = Interpolation.indminn(mod.mgrid.z, zmax, 1)[1]
	ixmin = Interpolation.indminn(mod.mgrid.x, xmin, 1)[1]
	ixmax = Interpolation.indminn(mod.mgrid.x, xmax, 1)[1]

	# calculate means
	mχvp = mean(mod.χvp);	mχvs = mean(mod.χvs);	mχρ = mean(mod.χρ)
	# remove means before Gaussian filtering
	mod.χvp -= mχvp; mod.χvs -= mχvs; mod.χρ -= mχρ
	Smooth.gaussian!(view(mod.χvp, izmin:izmax, ixmin:ixmax),[znwin,xnwin])
	Smooth.gaussian!(view(mod.χρ, izmin:izmax, ixmin:ixmax),[znwin,xnwin])
	Smooth.gaussian!(view(mod.χvs, izmin:izmax, ixmin:ixmax),[znwin,xnwin])
	# add mean back
	mod.χvp += mχvp; mod.χvs += mχvs; mod.χρ += mχρ
end

"""
Return a truncated seismic model using input bounds.
Note that there is no interpolation going on here, but only truncation, so 
the input bounds cannot be strictly imposed.

# Arguments
* `mod::Seismic` : model that is truncated

# Keyword Arguments

* `zmin::Float64=mod.mgrid.z[1]` : 
* `zmax::Float64=mod.mgrid.z[end]` : 
* `xmin::Float64=mod.mgrid.x[1]` : 
* `xmax::Float64=mod.mgrid.x[end]` : 
"""
function Seismic_trun(mod::Seismic;
			 zmin::Float64=mod.mgrid.z[1], zmax::Float64=mod.mgrid.z[end],
			 xmin::Float64=mod.mgrid.x[1], xmax::Float64=mod.mgrid.x[end] 
			 )

	izmin = Interpolation.indminn(mod.mgrid.z, zmin, 1)[1]
	izmax = Interpolation.indminn(mod.mgrid.z, zmax, 1)[1]
	ixmin = Interpolation.indminn(mod.mgrid.x, xmin, 1)[1]
	ixmax = Interpolation.indminn(mod.mgrid.x, xmax, 1)[1]

	x = mod.mgrid.x[ixmin:ixmax]
	z = mod.mgrid.z[izmin:izmax]
	npml =  mod.mgrid.npml
	δx, δz = mod.mgrid.δx, mod.mgrid.δz

	mgrid_trun = Grid.M2D(x, z, length(x), length(z), npml, δx, δz)
	
	# allocate model
	mod_trun = Seismic_zeros(mgrid_trun)

	for f in [:χvp, :χvs, :χρ]
		f0 = Symbol((replace("$(f)", "χ", "")),0)
		setfield!(mod_trun, f, getfield(mod, f)[izmin:izmax, ixmin:ixmax])
		setfield!(mod_trun, f0, getfield(mod, f0)) # adjust bounds later after truncation if necessary
	end
	return mod_trun

end

"""
Extend a seismic model into PML layers
"""
function Seismic_pml_pad_trun(mod::Seismic)
	modex=Seismic_zeros(Grid.M2D_pml_pad_trun(mod.mgrid))
	copy!(modex.vp0, mod.vp0)
	copy!(modex.vs0, mod.vs0)
	copy!(modex.ρ0, mod.ρ0)
	Seismic_pml_pad_trun!(modex, mod)

	return modex
end

"only padding implemented"
function Seismic_pml_pad_trun!(modex::Seismic, mod::Seismic)
	vp=mod.χvp; χ!(vp,mod.vp0,-1)
	vs=mod.χvs; χ!(vs,mod.vs0,-1)
	ρ=mod.χρ; χ!(ρ,mod.ρ0,-1)

	vpex=modex.χvp; vsex=modex.χvs; ρex=modex.χρ;
	pml_pad_trun!(vpex,vp,1);	pml_pad_trun!(vsex,vs,1);	pml_pad_trun!(ρex,ρ,1)
	χ!(vpex,modex.vp0,1);	χ!(vsex,modex.vs0,1);	χ!(ρex,modex.ρ0,1)
	χ!(vp,mod.vp0,1);	χ!(vs,mod.vs0,1);	χ!(ρ,mod.ρ0,1)
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
function pml_pad_trun!(modex::Array{Float64,2}, mod::Array{Float64,2}, flag::Int64=1)
	np=(size(modex,1)-size(mod,1) == size(modex,2)-size(mod,2)) ? 
			div(size(modex,2)-size(mod,2),2): error("modex size")
	if(isequal(flag,1)) 
		nz = size(mod,1); nx = size(mod,2)
		modexc=view(modex,np+1:np+nz,np+1:np+nx)
		@. modexc = mod
		for ix in 1:size(modex,2)
			for iz in 1:np
				modex[iz,ix] = modex[np+1,ix]
			end
			for iz in nz+1+np:size(modex,1)
				modex[iz,ix] = modex[nz+np,ix];
			end
		end
		for iz in 1:size(modex,1)
			for ix in 1:np
				modex[iz,ix] = modex[iz,np+1]
			end
			for ix in nx+1+np:size(modex,2)
				modex[iz,ix] = modex[iz,nx+np]
			end
		end
		return modex
	elseif(isequal(flag,-1)) 
		nz = size(modex,1); nx = size(modex,2)
		modexc=view(modex,np+1:nz-np,np+1:nx-np)
		@. mod=modexc
	else
		error("invalid flag")
	end
end


"""
function to resample in the model domain

# Arguments
* `mod::Seismic` : model
* `modi::Seismic` : model after interpolation
"""
function interp_spray!(mod::Seismic, mod_out::Seismic, attrib::Symbol, Battrib::Symbol=:B2; pa=nothing)
	if(pa===nothing)
		pa=Interpolation.Param([mod.mgrid.x, mod.mgrid.z], [mod_out.mgrid.x, mod_out.mgrid.z], Battrib)
	end

	"loop over fields in `Seismic`"
	Interpolation.interp_spray!(mod.χvp, mod_out.χvp, pa, attrib)
	Interpolation.interp_spray!(mod.χvs, mod_out.χvs, pa, attrib)
	Interpolation.interp_spray!(mod.χρ, mod_out.χρ,  pa, attrib)

	mod_out.vp0 = copy(mod.vp0); mod_out.vs0 = copy(mod.vs0); 
	mod_out.ρ0 = copy(mod.ρ0);
end
#
#macro inter

function save(mod::Seismic, folder)
	!(isdir(folder)) && error("invalid directory")

	nx=mod.mgrid.nx
	nz=mod.mgrid.nz
	n=max(nx,nz);
	fact=(n>40) ? round(Int,n/40) : 1
	mgrid=Grid.M2D_resamp(mod.mgrid, mod.mgrid.δx*fact, mod.mgrid.δz*fact)
	x=mgrid.x
	z=mgrid.z
	nx=length(x)
	nz=length(z)
	modo=Seismic_zeros(mgrid)
	interp_spray!(mod, modo, :interp)

	for m in [:vp, :ρ, :Zp]
		# save original gf
		file=joinpath(folder, string("im", m, ".csv"))
		CSV.write(file,DataFrame(hcat(repeat(z,outer=nx),
					      repeat(x,inner=nz),vec(Seismic_get(modo, m)))))
	end
end
	

end # module
