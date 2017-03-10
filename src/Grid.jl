module Grid

type M2D
	x::Array{Float64}
	z::Array{Float64}
	nx::Int64
	nz::Int64
	npml::Int64
	δx::Float64
	δz::Float64
end

function(attrib::Symbol)
	if(attrib == :samp1)
		return M2D(-1000.0,1000.0,-1000.0,1000.0,201,201,40)
	else
		error("invalid attrib")
	end
end

function M2D(xmin::Float64, xmax::Float64,
	zmin::Float64, zmax::Float64,
	nx::Int64, nz::Int64,
	npml::Int64 
	)
	x = Array(linspace(xmin,xmax,nx));
	z = Array(linspace(zmin,zmax,nz));
	return M2D(x, z, nx, nz, npml, x[2]-x[1], z[2]-z[1])
end

function M2D(xmin::Float64, xmax::Float64,
	zmin::Float64, zmax::Float64,
	δx::Float64, δz::Float64,
	npml::Int64,
	attrib::Symbol=nothing
	)
	x = [xx for xx in xmin:δx:xmax]
	z = [zz for zz in zmin:δz:zmax]
	nx = size(x,1); nz = size(x,1);
	return M2D(x, z, nx, nz, npml, x[2]-x[1], z[2]-z[1])
end


type M1D
	x::Array{Float64}
	nx::Int64
	δx::Float64
#	OrderedPair(x,y) = x > y ? error("out of order") : new(x,y)
end

"""
Default M1D grids based on attirb
"""
function M1D(attrib::Symbol)
	if(attrib == :timesamp1)
		return M1D(0.0,1.0,2000)
	elseif(attrib == :npow2samp)
		return M1D(npow2=16,δ=0.0001)
	else
		error("invalid attrib")
	end
end

"""
Construct 1-D grid based on nx or δx
"""
function M1D(xbeg::Float64, xend::Float64, nx::Int64)
	x = Array(linspace(xbeg, xend, nx))
	return M1D(x, nx, x[2]-x[1])
end

function M1D(xbeg::Float64, xend::Float64, δx::Float64)
	x = [tt for tt in xbeg:δx:xend]
	return M1D(x, size(x,1), x[2]-x[1])
end


"""
output an npow2grid of either time or frequency
"""
function M1D_npow2(npow2::Int64, δ::Float64)
	vec = zeros(npow2);

	# zero lag
	vec[1] = 0.0;

	# positive lags
	for i = 1: div(npow2,2)
		vec[1+i] = δ * i
	end

	# negative lags -- one less than positive lags
	for i = 1: div(npow2,2)-1
		vec[npow2-i+1] = -δ * i
	end
	return M1D(vec, npow2, δ)
end

"convertion between time and frequency npow2 grids"
M1D_npow2_tf(grid::M1D) = M1D_npow2(grid.nx, 1.0/grid.nx/grid.δx)



end # module
