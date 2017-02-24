module Grid

type Mod2D
	x::Array{Float64}
	z::Array{Float64}
	nx::Int64
	nz::Int64
	npml::Int64
	δx::Float64
	δz::Float64
end

function Mod2D()
	nx = 200; nz = 200; npml = 40;
	x = linspace(-1000,1000.0,nx);
	z = linspace(-1000,1000.0,nz);
	return Mod2D(x, z, nx, nz, npml, x[2]-x[1], z[2]-z[1])
end

type Tim
	t::Array{Float64}
	nt::Int64
	δt::Float64
end

function Tim()
	nt = 2000;
	t = linspace(0.0,1.0,2000)
	return Tim(t, nt, t[2]-t[1])
end

end # module
