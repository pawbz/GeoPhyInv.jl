module Models

import SIT.Grid: M2D

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
	mgrid::M2D
	attrib::Symbol
end

function Seismic(attrib::Symbol)
	if(attrib == :homo_acoustic1)
		vp0 = 2000.0;
		vs0 = 0.0;
		ρ0 = 2000.0;
		mgrid = M2D(:samp1)
		return Seismic(vp0, vs0, ρ0,
		      fill(vp0, (mgrid.nz, mgrid.nx)),
		      fill(vs0, (mgrid.nz, mgrid.nx)),
		      fill(ρ0, (mgrid.nz, mgrid.nx)),
		      mgrid, attrib)
	end
end

function Seismic(vp0::Float64,
		 vs0::Float64,
		 ρ0::Float64,
		 vp::Array{Float64},
		 vs::Array{Float64},
		 ρ::Array{Float64},
		 mgrid::M2D,
		 attrib::AbstractString)

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
	      mgrid,
	      attrib)
end

function χ(mod::Array{Float64}, mod0::Float64, flag::Int64=1)
	if(flag == 1)
		return	((mod - mod0) * mod0^(-1.0))
	elseif(flag == -1)
		return  (mod .* mod0 + mod0)
	end
end


end # module
