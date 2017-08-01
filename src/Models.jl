__precompile__()

module Models


import JuMIT.Grid
import JuMIT.IO
import JuMIT.Interpolation
import JuMIT.Smooth

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
function print(mod::Seismic, name::String="")
	println("\tSeismic Model:\t",name)
	println("\t> number of samples:\t","x\t",mod.mgrid.nx,"\tz\t",mod.mgrid.nz)
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

function adjust_bounds!(mod::Seismic, frac::Float64)
	for f in [:χvp, :χρ, :χvs]
		f0 = Symbol((replace("$(f)", "χ", "")),0)
		fm = Symbol((replace("$(f)", "χ", "")))
		m = Seismic_get(mod, fm)
		setfield!(mod, f0, bounds(m , frac))
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

function Seismic_iszero(mod::Seismic)
	return any((mod.vp0 .* mod.ρ0) .== 0.0) ? true : false
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
	K0 = mod.vp0 .* mod.vp0 .* mod.ρ0; 
	KI0 = flipdim(K0.^(-1.0),1);
	"note that: mean(K0) ≠ 1/mean(KI0)"
	ρ0 = mod.ρ0; ρI0 = flipdim(ρ0.^(-1.0),1);
	μ0 = mod.vs0 .* mod.vs0 .* mod.ρ0;
	vp = χ(mod.χvp, mod.vp0, -1)
	vs = χ(mod.χvs, mod.vs0, -1)
	ρ = χ(mod.χρ, mod.ρ0, -1)
	if(attrib == :ρI)
		return ρ.^(-1.0);
	elseif(attrib == :ρ)
		return ρ;
	elseif(attrib == :vp)
		return vp;
	elseif(attrib == :vs)
		return vs;
	elseif(attrib == :ρI0)
		return ρI0
	elseif(attrib == :χρI)
		return χ(ρ.^(-1.0), ρI0);
	elseif(attrib == :χK)
		return χ(vp .* vp .* ρ, K0);
	elseif(attrib == :χμ)
		return χ(vs .* vs .* ρ, μ0);
	elseif(attrib == :K0)
		return K0
	elseif(attrib == :μ0)
		return μ0
	elseif(attrib == :KI0)
		return KI0
	elseif(attrib == :χKI)
		return χ((vp .* vp .* ρ).^(-1), KI0); 
	elseif(attrib == :KI)
		return (vp .* vp .* ρ).^(-1);
	elseif(attrib == :K)
		return (vp .* vp .* ρ);
	elseif(attrib == :Zp)
		return (vp .* ρ);
	else
		error("invalid attrib")
	end
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
	Seismic_iszero(mod) ? error("mod cannot be zero") : nothing
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
		      g1::Vector{Float64},
		      g2::Vector{Float64},
		      attribs::Vector{Symbol}=[:χKI, :χρI],
		      flag::Int64=1
		      )

	if(flag == 1)
		if(attribs == [:χKI, :χρI]) 
			vp0 = mod.vp0;	vs0 = mod.vs0;	ρ0 = mod.ρ0
			ρI = Seismic_get(mod,:ρI);	vp = Seismic_get(mod,:vp)
			Zp = Seismic_get(mod,:Zp)

			gKI = χg(reshape(g1, mod.mgrid.nz, mod.mgrid.nx), Seismic_get(mod,:KI0), -1)
			gρI = χg(reshape(g2, mod.mgrid.nz, mod.mgrid.nx), Seismic_get(mod,:ρI0), -1)

			# gradient w.r.t  χρ
			gmod.χρ = copy(χg((-1. .* ρI.^2 .* gρI - (Zp).^(-2) .* gKI), ρ0, 1))
			# gradient w.r.t χvp
			gmod.χvp = copy(χg((-2. .* (vp).^(-3) .* ρI .* gKI), vp0, 1))
			# to be implemented later
			gmod.χvs = copy(zeros(gmod.χvp))
			gmod.mgrid = deepcopy(mod.mgrid);	gmod.vp0 = copy(mod.vp0);	gmod.ρ0 = copy(mod.ρ0)
		else
			error("invalid attribs")
		end
	elseif(flag == -1)
		gvp = vec(χg(gmod.χvp, mod.vp0, -1))
		gρ = vec(χg(gmod.χρ, mod.ρ0, -1))
		K0 = mean(Seismic_get(mod, :K0))
		KI0 = mean(Seismic_get(mod, :KI0))
		ρ0 = mean(mod.ρ0)
		ρI0 = mean(Seismic_get(mod, :ρI0))
		ρI = vec(Seismic_get(mod,:ρI));	vp = vec(Seismic_get(mod,:vp))
		ρ = vec(Seismic_get(mod,:ρ)); ρ0 = mean(mod.ρ0)

		nznx = mod.mgrid.nz * mod.mgrid.nx;
		if(attribs == [:χK, :χρ])
			g1[1:nznx] = copy(0.5 .* K0 * gvp .* ρI .* vp.^(-1))
			g2[1:nznx] = copy((gρ .* ρ0) - (0.5 .* ρ0 .* (gvp) .* (vp) ./ (ρ) ))
		elseif(attribs == [:χKI, :χρI])
			g1[1:nznx] =  copy(-0.5 .* KI0 .* gvp .* ρ .*  vp.^(3.))
			g2[1:nznx] = copy((gρ .* ρI0 .* (-1.) .*  (ρ.^(2))) +
			    (0.5 * ρI0 .* (gvp) .*  (vp) .* (ρ)))
		else
			error("invalid attribs")
		end
	else
		error("invalid flag")
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
function χ(mod::Array{Float64}, mod0::Vector{Float64}, flag::Int64=1)
	m0 = mean(mod0)
	if(flag == 1)
		# return zero when m0 is zero
		return	(m0 == 0.0) ? 0.0 : ((mod - m0) * m0^(-1.0))
	elseif(flag == -1)
		return  (mod .* m0 + m0)
	end
end

"""
Gradients
Return contrast model parameter using the
reference value.
"""
function χg(mod::Array{Float64}, mod0::Vector{Float64}, flag::Int64=1)
	m0 = mean(mod0)
	if(flag == 1)
		return  mod * m0
	elseif(flag == -1)
		# return zero when m0 is zero
		return	(m0 == 0.0) ? 0.0 : mod * m0^(-1.0)
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
* `randn_pert::Float64=0.0` : percentage of reference values for additive random noise
"""
function Seismic_addon!(mod::Seismic; 
		       circ_loc::Vector{Float64}=[0., 0.,],
		       circ_rad::Float64=0.0,
		       circ_pert::Float64=0.0,
		       rect_loc::Vector{Float64}=[0., 0., 0., 0.,],
		       rect_pert::Float64=0.0,
		       randn_perc::Real=0.0,
		       fields::Vector{Symbol}=[:χvp,:χρ,:χvs]
		       )

	temp = zeros(mod.mgrid.nz, mod.mgrid.nx)
	if(!(circ_pert == 0.0))
		temp += [((sqrt((mod.mgrid.x[ix]-circ_loc[2])^2 + 
			(mod.mgrid.z[iz]-circ_loc[1])^2 ) <= circ_rad) ? circ_pert : 0.0)  for 
	   		iz in 1:mod.mgrid.nz, ix in 1:mod.mgrid.nx ]
	end
	if(!(rect_pert == 0.0))
		temp += [
			(((mod.mgrid.x[ix]-rect_loc[4]) * (mod.mgrid.x[ix]-rect_loc[2]) < 0.0) & 
			((mod.mgrid.z[iz]-rect_loc[3]) * (mod.mgrid.z[iz]-rect_loc[1]) < 0.0)) ?
			rect_pert : 0.0  for
			iz in 1:mod.mgrid.nz, ix in 1:mod.mgrid.nx ]
	end

	for field in fields
		setfield!(mod, field, (getfield(mod, field)+temp))
	end
	# random noise (in future change to fields)
	mod.χvp += randn(size(mod.χvp)) .* randn_perc .* 1e-2 .* mean(mod.vp0); 
	mod.χvs += randn(size(mod.χvs)) .* randn_perc .* 1e-2 .* mean(mod.vs0); 
	mod.χρ += randn(size(mod.χρ)) .* randn_perc .* 1e-2 .* mean(mod.ρ0); 

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
PML Extend a model on all four sides
"""
function pad_trun(mod::Array{Float64,2}, np::Int64, flag::Int64=1)
	if(isequal(flag,1)) 
		nz = size(mod,1); nx = size(mod,2)
		modex = zeros(nz + 2*np, nx + 2*np)

		modex[np+1:np+nz,np+1:np+nx] .= mod 
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
function Seismic_interp_spray!(mod::Seismic, mod_out::Seismic, attrib::Symbol, Battrib::Symbol=:B2 )

	"loop over fields in `Seismic`"
	Interpolation.interp_spray!(mod.mgrid.x, mod.mgrid.z, mod.χvp,
		      mod_out.mgrid.x, mod_out.mgrid.z, mod_out.χvp, attrib, Battrib)
	Interpolation.interp_spray!(mod.mgrid.x, mod.mgrid.z, mod.χvs,
		      mod_out.mgrid.x, mod_out.mgrid.z, mod_out.χvs, attrib, Battrib)
	Interpolation.interp_spray!(mod.mgrid.x, mod.mgrid.z, mod.χρ,
		      mod_out.mgrid.x, mod_out.mgrid.z, mod_out.χρ, attrib, Battrib)

	mod_out.vp0 = copy(mod.vp0); mod_out.vs0 = copy(mod.vs0); 
	mod_out.ρ0 = copy(mod.ρ0);
end
#
#macro inter


end # module
