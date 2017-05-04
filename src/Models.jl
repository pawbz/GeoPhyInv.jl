module Models


import SIT.Grid
import SIT.IO
using Interpolations

"""
Data type fo represent a seismic model.
# Fields
* `vp0` :
* `vs0` :
* `ρ0` :
* `χvp` :
* `χvs` :
* `χρ` :
* `mgrid` :
"""
type Seismic
	vp0::Float64
	vs0::Float64
	ρ0::Float64
	χvp::Array{Float64}
	χvs::Array{Float64}
	χρ::Array{Float64}
	mgrid::Grid.M2D
	"adding conditions that are to be false while construction"
	Seismic(vp0,vs0,ρ0,χvp,χvs,χρ,mgrid) = 
		any([
       	 	     vp0<0.0, 
		     vs0<0.0,
		     ρ0<0.0,
		     size(χvp) != (length(mgrid.z), length(mgrid.x)),
		     size(χvs) != (length(mgrid.z), length(mgrid.x)),
		     size(χρ) != (length(mgrid.z), length(mgrid.x))
		    ]) ? 
		       error("Seismic construct") : new(vp0,vs0,ρ0,χvp,χvs,χρ,mgrid)
end

"""
Return `Seismic` with zeros everywhere;
this method is used for preallocation.

# Arguments
* `mgrid::Grid.M2D` : 
"""
function Seismic_zeros(mgrid::Grid.M2D)
	return Seismic(0.0, 0.0, 0.0, zeros(mgrid.nz, mgrid.nx), zeros(mgrid.nz, mgrid.nx),
		zeros(mgrid.nz, mgrid.nx), mgrid)
end

function Seismic_iszero(mod::Seismic)
	return mod.vp0 * mod.ρ0 == 0.0 ? true : false
end

"Logical operation for `Seismic`"
function Seismic_isequal(mod1::Seismic, mod2::Seismic)
	fnames = fieldnames(Seismic)
	pop!(fnames) # last one has different isequal
	vec=[(isequal(getfield(mod1, name),getfield(mod2, name))) for name in fnames]
	push!(vec, Grid.M2D_isequal(mod1.mgrid, mod2.mgrid))
	return all(vec)
end

"Logical operation for `Seismic`"
function Seismic_issimilar(mod1::Seismic, mod2::Seismic)
	return Grid.M2D_isequal(mod1.mgrid, mod2.mgrid)
end

"Logical operation for `Seismic`"
function Seismic_issimilar(mod1::Seismic, grid::Grid.M2D)
	return Grid.M2D_isequal(mod1.mgrid, grid)
end



"""
Get other dependent model parameters of a seismic model
that are not present in `Seismic`.

* `:ρI` : inverse of density
* `:Zp` : P-wave impedance
"""
function Seismic_get(mod::Seismic, attrib::Symbol)
	K0 = mod.vp0 * mod.vp0 * mod.ρ0; K0I = K0^(-1.0);
	μ0 = mod.vs0 * mod.vs0 * mod.ρ0;
	vp = χ(mod.χvp, mod.vp0, -1)
	vs = χ(mod.χvs, mod.vs0, -1)
	ρ = χ(mod.χρ, mod.ρ0, -1)
	if(attrib == :ρI)
		return ρ.^(-1.0)
	elseif(attrib == :ρ)
		return ρ
	elseif(attrib == :vp)
		return vp
	elseif(attrib == :ρ0I)
		return mod.ρ0.^(-1.0)
	elseif(attrib == :χρI)
		return χ(ρ.^(-1.0), mod.ρ0^(-1.0))
	elseif(attrib == :χK)
		return χ(vp .* vp .* ρ, K0) 
	elseif(attrib == :χμ)
		return χ(vs .* vs .* ρ, μ0) 
	elseif(attrib == :K0)
		return K0
	elseif(attrib == :μ0)
		return μ0
	elseif(attrib == :K0I)
		return K0I
	elseif(attrib == :χKI)
		return χ((vp .* vp .* ρ).^(-1), K0I) 
	elseif(attrib == :KI)
		return (vp .* vp .* ρ).^(-1)
	elseif(attrib == :K)
		return (vp .* vp .* ρ)
	elseif(attrib == :Zp)
		return (vp .* ρ)
	else
		error("invalid attrib")
	end
end

"""
Use chain rule to output gradients after
re-parameterization.

# Arguments
* `gmod::Seismic` : gradient model 
* `mod::Seismic` : model 
* `attribs::Vector{Symbol}` : parameterization, for example [:χKI, :χρI]
* `g1::Vector{Float64}` : 
* `g2::Vector{Float64}` : 
"""
function Seismic_getgrad!(gmod::Seismic, 
			  mod::Seismic, 
			  attribs::Vector{Symbol}, 
			  g1::Vector{Float64}, g2::Vector{Float64})
	gvp = vec(χg(gmod.χvp, mod.vp0, -1))
	gρ = vec(χg(gmod.χρ, mod.ρ0, -1))
	K0 = Seismic_get(mod, :K0)
	ρI = vec(Seismic_get(mod,:ρI));	vp = vec(Seismic_get(mod,:vp))
	ρ = vec(Seismic_get(mod,:ρ)); ρ0 = mod.ρ0 

	nznx = mod.mgrid.nz * mod.mgrid.nx;
	if(attribs == [:χK, :χρ])
		g1[1:nznx] = copy(0.5 .* K0 * gvp .* ρI .* vp.^(-1))
		g2[1:nznx] = copy((gρ .* ρ0) - (0.5 .* ρ0 .* (gvp) .* (vp) ./ (ρ) ))
	elseif(attribs == [:χKI, :χρI])
		g1[1:nznx] =  copy(-0.5 .* K0^(-1) .* gvp .* ρ .*  vp.^(3.))
		g2[1:nznx] = copy((gρ .* ρ0^(-1) .* (-1.) .*  (ρ.^(2))) +
		    (0.5 * ρ0^(-1.) .* (gvp) .*  (vp) .* (ρ)))
	else
		error("invalid attribs")
	end
	return g1, g2
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
	if(attribs = [:χKI, :χρI]) 
		ρ = (χ(x2, Seismic_get(mod, :ρ0I), -1)).^(-1)
		K = (χ(x1, Seismic_get(mod, :K0I), -1)).^(-1)
		vp = broadcast(sqrt, (K .* (χ(x2, Seismic_get(mod, :ρ0I), -1))))
		mod.χvp = copy(χ(vp, mod.vp0, 1))
		mod.χρ = copy(χ(ρ, mod.ρ0, 1))
	else
		error("invalid attribs")
	end
end

"""
Use chain rule to output gradients with 
respect to χvp and χρ from  gradients 
with respect to KI and ρI.

# Arguments
* `mod::Seismic` : model required for chain rule
* `gKI` : gradient of an objective function with respect to KI
* `gρI` : gradient of an objective function with respect to ρI
"""
function Seismic_chainrule!(
		      gmod::Seismic,
		      mod::Seismic,
		      g1::Vector{Float64},
		      g2::Vector{Float64},
		      attribs::Vector{Symbol}=[:χKI, :χρI]
		      )

	if(attribs = [:χKI, :χρI]) 
		vp0 = mod.vp0;	vs0 = mod.vs0;	ρ0 = mod.ρ0
		ρI = Seismic_get(mod,:ρI);	vp = Seismic_get(mod,:vp)
		Zp = Seismic_get(mod,:Zp)

		gKI = χg(reshape(g1, mod.mgrid.nz, mod.mgrid.nx), Seismic_get(mod,:K0I), -1)
		gρI = χg(reshape(g2, mod.mgrid.nz, mod.mgrid.nx), Seismic_get(mod,:ρ0I), -1)

		# gradient w.r.t  χρ
		gmod.χρ = copy(χg((-1. .* ρI.^2 .* gρI - (Zp).^(-2) .* gKI), ρ0, 1))
		# gradient w.r.t χvp
		gmod.χvp = copy(χg((-2. .* (vp).^(-3) .* ρI .* gKI), vp0, 1))
		# to be implemented later
		gmod.χvs = copy(zeros(gmod.χvp))
		gmod.mgrid = deepcopy(mod.mgrid);	gmod.vp0 = copy(mod.vp0);	gmod.ρ0 = copy(mod.ρ0)

		return gmod
	else
		error("invalid attribs")
	end
end


"""
Return dimensionless contrast model parameter using the
reference value.

# Arguments
* `mod::Array{Float64}` : subsurface parameter
* `mod0::Float64` : reference value
* `flag::Int64=1` : 
"""
function χ(mod::Array{Float64}, mod0::Float64, flag::Int64=1)
	if(flag == 1)
		return	((mod - mod0) * mod0^(-1.0))
	elseif(flag == -1)
		return  (mod .* mod0 + mod0)
	end
end

"""
Gradients
Return contrast model parameter using the
reference value.
"""
function χg(mod::Array{Float64}, mod0::Float64, flag::Int64=1)
	if(flag == 1)
		return  mod * mod0
	elseif(flag == -1)
		return	mod * mod0^(-1.0)
	end
end


"""
Add features to a model.

# Arguments
* `mod::Seismic` : model that is modified

# Keyword Arguments
* `circ_loc::Vector{Float64}=nothing` : location of center of perturbation
* `circ_rad::Float64=0.0` : radius of circular perturbation
* `circ_pert::Float64=0.1` : perturbation inside a circle
* `rect_loc::Array{Float64}=nothing` : rectangle location
* `rect_pert::Float64=0.1` : perturbation in a rectangle

"""
function Seismic_addon(mod::Seismic; 
		       circ_loc::Vector{Float64}=[0., 0.,],
		       circ_rad::Float64=0.0,
		       circ_pert::Float64=0.1,
		       rect_loc::Vector{Float64}=[0., 0., 0., 0.,],
		       rect_pert::Float64=0.1,
		       fields::Vector{Symbol}=[:χvp,:χρ,:χvs]
		       )

	temp = zeros(mod.mgrid.nz, mod.mgrid.nx)
	if(!(circ_loc == nothing))
		temp += [((sqrt((mod.mgrid.x[ix]-circ_loc[2])^2 + 
			(mod.mgrid.z[iz]-circ_loc[1])^2 ) <= circ_rad) ? circ_pert : 0.0)  for 
	   		iz in 1:mod.mgrid.nz, ix in 1:mod.mgrid.nx ]
	end
	if(!(rect_loc == nothing))
		temp += [
			(((mod.mgrid.x[ix]-rect_loc[4]) * (mod.mgrid.x[ix]-rect_loc[2]) < 0.0) & 
			((mod.mgrid.z[iz]-rect_loc[3]) * (mod.mgrid.z[iz]-rect_loc[1]) < 0.0)) ?
			rect_pert : 0.0  for
			iz in 1:mod.mgrid.nz, ix in 1:mod.mgrid.nx ]
	end

	mod_out = deepcopy(mod)
	for field in fields
		setfield!(mod_out, field, (getfield(mod, field)+temp))
	end
	return mod_out
end

"""
Extend a seismic model into PML layers
"""
function Seismic_pad_trun(mod::Seismic)

	vpex = pad_trun(χ(mod.χvp, mod.vp0, -1),mod.mgrid.npml);
	vsex = pad_trun(χ(mod.χvs, mod.vs0, -1),mod.mgrid.npml);
	ρex = pad_trun(χ(mod.χρ, mod.ρ0, -1),mod.mgrid.npml);
	return Seismic(mod.vp0,
		 mod.vs0,
		 mod.ρ0,
		 χ(vpex,mod.vp0),
		 χ(vsex,mod.vs0),
		 χ(ρex,mod.ρ0),
		 Grid.M2D_pad_trun(mod.mgrid))
end

"""
Extend a model on all four sides
"""
function pad_trun(mod::Array{Float64,2}, np::Int64, flag::Int64=1)
	if(isequal(flag,1)) 
		nz = size(mod,1); nx = size(mod,2)
		modex = zeros(nz + 2*np, nx + 2*np)

		modex[np+1:np+nz,np+1:np+nx] = mod 
		modex[1:np,:] = repeat(modex[np+1,:], inner=np);
		modex[nz+1+np:end,:] = repeat(modex[nz+np,:],inner=np);
		modex[:,1:np] = repeat(modex[:,np+1], outer=np);
		modex[:,nx+1+np:end] = repeat(modex[:,nx+np], outer=np);
		return modex
	elseif(isequal(flag,-1)) 
		nz = size(mod,1); nx = size(mod,2)
		return mod[np+1:nz-np,np+1:nx-np]
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
function Seismic_interp_spray!(mod::Seismic, mod_out::Seismic, attrib::Symbol)

	itpvp = interpolate((mod.mgrid.z, mod.mgrid.x),
		     mod.χvp, 
		     Gridded(Linear()))
	itpvs = interpolate((mod.mgrid.z, mod.mgrid.x),
		     mod.χvs, 
		     Gridded(Linear()))

	itpρ = interpolate((mod.mgrid.z, mod.mgrid.x),
		     mod.χρ, 
		     Gridded(Linear()))

	if(attrib == :interp) 
		mod_out.χvp = copy(itpvp[mod_out.mgrid.z, mod_out.mgrid.x])
		mod_out.χρ = copy(itpρ[mod_out.mgrid.z, mod_out.mgrid.x])
		mod_out.χvs = copy(itpvs[mod_out.mgrid.z, mod_out.mgrid.x])

	elseif(attrib == :spray)

		mod_out.χvp = copy(gradient(itpvp[mod_out.mgrid.z, mod_out.mgrid.x]))
		mod_out.χρ = copy(gradient(itpρ[mod_out.mgrid.z, mod_out.mgrid.x]))
		mod_out.χvs = copy(gradient(itpvs[mod_out.mgrid.z, mod_out.mgrid.x]))

	else
		error("invalid flag")
	end
	mod_out.vp0 = copy(mod.vp0); mod_out.vs0 = copy(mod.vs0); 
	mod_out.ρ0 = copy(mod.ρ0);
end
#
#macro inter


end # module
