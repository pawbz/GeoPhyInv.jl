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
* `ρ0I` :
* `K0` :
* `K0I` :
* `μ0` :
* `χvp` :
* `χvs` :
* `χρ` :
* `χρI` :
* `χK` :
* `χKI` :
* `χμ` :
* `mgrid` :
"""
type Seismic
	vp0::Float64
	vs0::Float64
	ρ0::Float64
	ρ0I::Float64
	K0::Float64
	K0I::Float64
	μ0::Float64
	χvp::Array{Float64}
	χvs::Array{Float64}
	χρ::Array{Float64}
	χρI::Array{Float64}
	χK::Array{Float64}
	χKI::Array{Float64}
	χμ::Array{Float64}
	mgrid::Grid.M2D
end

"""
Construct using input velocity and density arrays
"""
function Seismic(vp0::Float64,
		 vs0::Float64,
		 ρ0::Float64,
		 vp::Array{Float64},
		 vs::Array{Float64},
		 ρ::Array{Float64},
		 mgrid::Grid.M2D)

	ρ0I = ρ0^(-1.0);
	K0 = vp0 * vp0 * ρ0; K0I = K0^(-1.0);
	μ0 = vs0 * vs0 * ρ0;
	return Seismic(vp0,vs0,ρ0,ρ0I,K0,K0I,μ0,
	      χ(vp, vp0),
	      χ(vs, vs0),
	      χ(ρ, ρ0),
	      χ(ρ.^(-1.0), ρ0I),
	      χ(vp .* vp .* ρ, K0),
	      χ((vp .* vp .* ρ).^(-1), K0I),
	      χ(vs .* vs .* ρ, μ0),
	      mgrid)
end


"""
Return contrast model parameter using the
reference value.
"""
function χ(mod::Array{Float64}, mod0::Float64, flag::Int64=1)
	if(flag == 1)
		return	((mod - mod0) * mod0^(-1.0))
	elseif(flag == -1)
		return  (mod .* mod0 + mod0)
	end
end


"""
Extend a seismic model into PML layers
"""
function Seismic_extend(mod::Seismic)

	vpex = extend(χ(mod.χvp, mod.vp0, -1),mod.mgrid.npml);
	vsex = extend(χ(mod.χvs, mod.vs0, -1),mod.mgrid.npml);
	ρex = extend(χ(mod.χρ, mod.ρ0, -1),mod.mgrid.npml);
	return Seismic(mod.vp0,
		 mod.vs0,
		 mod.ρ0,
		 vpex,
		 vsex,
		 ρex,
		 Grid.M2D_extend(mod.mgrid))
end

"""
Extend a model on all four sides
"""
function extend(mod::Array{Float64,2}, np::Int64)
	nz = size(mod,1); nx = size(mod,2)
	modex = zeros(nz + 2*np, nx + 2*np)

	modex[np+1:np+nz,np+1:np+nx] = mod 
	modex[1:np,:] = repeat(modex[np+1,:], inner=np);
	modex[nz+1+np:end,:] = repeat(modex[nz+np,:],inner=np);
	modex[:,1:np] = repeat(modex[:,np+1], outer=np);
	modex[:,nx+1+np:end] = repeat(modex[:,nx+np], outer=np);
	return modex
end


"""
function to resample in the model domain
"""
function Seismic_resamp(mod::Seismic,
		mgrid::Grid.M2D
		)

	itpvp = interpolate((mod.mgrid.z, mod.mgrid.x),
		     mod.χvp, 
		     Gridded(Linear()))
	itpvs = interpolate((mod.mgrid.z, mod.mgrid.x),
		     mod.χvs, 
		     Gridded(Linear()))

	itpρ = interpolate((mod.mgrid.z, mod.mgrid.x),
		     mod.χρ, 
		     Gridded(Linear()))

	return Seismic(mod.vp0, mod.vs0, mod.ρ0, 
		χ(itpvp[mgrid.z, mgrid.x], mod.vp0, -1), 
		χ(itpvs[mgrid.z, mgrid.x], mod.vs0, -1),
		χ(itpρ[mgrid.z, mgrid.x], mod.ρ0, -1),
			      mgrid)

end
#
#macro inter


end # module
