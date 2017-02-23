module models

export SeismicModel, χ

using SeismicInversion.mesh

type SeismicModel
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
	attrib::AbstractString
end

function SeismicModel(attrib::AbstractString, mesh::Mesh2D)
	if(attrib == "test_homo_acoustic")
		vp0 = 2000.0;
		vs0 = 0.0;
		ρ0 = 2000.0; 
		return SeismicModel(vp0, vs0, ρ0,
		      vp0 .* ones(mesh.nz, mesh.nx),
		      vs0 .* ones(mesh.nz, mesh.nx),
		      ρ0 .* ones(mesh.nz, mesh.nx),
		      attrib, mesh)
	end
end

function SeismicModel(vp0::Float64, 
		      vs0::Float64, 
		      ρ0::Float64, 
		      vp::Array{Float64}, 
		      vs::Array{Float64}, 
		      ρ::Array{Float64},
		      attrib::AbstractString, mesh::Mesh2D)
		ρ0I = ρ0^(-1.0);
		K0 = vp0 * vp0 * ρ0; K0I = K0^(-1.0);
		μ0 = vs0 * vs0 * ρ0;
		return SeismicModel(vp0,vs0,ρ0,ρ0I,K0,K0I,μ0,
		      χ(vp, vp0),
		      χ(vs, vs0),
		      χ(ρ, ρ0),
		      χ(ρ.^(-1.0), ρ0I),
		      χ(vp .* vp .* ρ, K0),
		      χ((vp .* vp .* ρ).^(-1), K0I),
		      χ(vs .* vs .* ρ, μ0),
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
