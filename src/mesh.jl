module mesh

export Mesh2D

type Mesh2D
	x :: Array{Float64}
	z :: Array{Float64}
	nx :: Int64
	nz :: Int64
	npml :: Int64
	δx :: Float64
	δz :: Float64
end

function Mesh2D()
	nx = 200; nz = 200; npml = 40;
	x = linspace(0.0,2000.0,nx);
	z = linspace(0.0,2000.0,nz);
	return Mesh2D(x, z, nx, nz, npml, x[2]-x[1], z[2]-z[1])
end

end # module
