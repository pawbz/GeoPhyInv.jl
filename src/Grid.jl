module Grid

"""
Data type to represent 2D grid.
# Fields
* `x` : horizontal values
* `z` : vertical values
* `nx` : number of samples in horizontal direction
* `nz` : number of samples in vertical direction
* `δx` : sampling interval in horizontal direction
* `δz` : sampling interval in vertical direction
"""

type M2D
	x::Array{Float64}
	z::Array{Float64}
	nx::Int64
	nz::Int64
	npml::Int64
	δx::Float64
	δz::Float64
	"adding conditions that are to be false while construction"
	M2D(x, z, nx, nz, npml, δx, δz) = 
		any([
       		  δx < 0.0, length(x) != nx,
       		  δz < 0.0, length(z) != nz
		  ]) ? 
		error("M2D construct") : new(x, z, nx, nz, npml, δx, δz)
end

"Logical operation for `M2D`"
function M2D_isequal(grid1::M2D, grid2::M2D)
	return all([(isequal(getfield(grid1, name),getfield(grid2, name))) for name in fieldnames(M2D)])
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
	)
	x = [xx for xx in xmin:δx:xmax]
	z = [zz for zz in zmin:δz:zmax]
	nx = size(x,1); nz = size(x,1);
	return M2D(x, z, nx, nz, npml, x[2]-x[1], z[2]-z[1])
end

"""
2-D grid with a different sampling interval
"""
M2D_resamp(grid::M2D, δx::Float64, δz::Float64) = M2D(grid.x[1], grid.x[end], 
					grid.z[1], grid.z[end], δx, δz, grid.npml)

"""
Extend M2D by its npml grid points on all sides
"""
function M2D_pad_trun(mgrid::M2D; flag::Int64=1)

	if(isequal(flag,1)) 
		xmin = mgrid.x[1] - mgrid.npml*mgrid.δx
		xmax = mgrid.x[end] + mgrid.npml*mgrid.δx
		zmin = mgrid.z[1] - mgrid.npml*mgrid.δz
		zmax = mgrid.z[end] + mgrid.npml*mgrid.δz
		return M2D(xmin,xmax,zmin,zmax, mgrid.nx+2*mgrid.npml,mgrid.nz+2*mgrid.npml,mgrid.npml)
	elseif(isequal(flag,-1))
		xmin = mgrid.x[1] + mgrid.npml*mgrid.δx
		xmax = mgrid.x[end] - mgrid.npml*mgrid.δx
		zmin = mgrid.z[1] + mgrid.npml*mgrid.δz
		zmax = mgrid.z[end] - mgrid.npml*mgrid.δz
		return M2D(xmin,xmax,zmin,zmax, mgrid.nx-2*mgrid.npml,mgrid.nz-2*mgrid.npml,mgrid.npml)
	else
		error("invalid flag")
	end
	
end

"""
Return the X and Z positions of the boundary of mgrid
attrib
** :inner
** :outer
** onlycount
"""
function M2D_boundary(mgrid::M2D, nlayer::Int64, attrib::Symbol; onlycount::Bool=false)
	if(attrib == :inner)
		x = mgrid.x; z = mgrid.z;
	elseif(attrib == :outer)
		# extending x and z
		x = vcat(
	   		[mgrid.x[1]-(ilayer)*mgrid.δx for ilayer=1:nlayer:-1],
	   		mgrid.x,
			[mgrid.x[end]+(ilayer)*mgrid.δx for ilayer=1:nlayer],
			)
		z = vcat(
	   		[mgrid.z[1]-(ilayer)*mgrid.δz for ilayer=1:nlayer:-1],
	   		mgrid.z,
			[mgrid.z[end]+(ilayer)*mgrid.δz for ilayer=1:nlayer],
			)
	end
	nx = length(x); nz = length(z);
	bx = vcat(
		  x, x[end]*ones(nz-1),	  x[end-1:-1:1], x[1]*ones(nz-2))
	bz = vcat(
		  z[1]*ones(nx), z[2:end],  z[end]*ones(nx-1), z[end-1:-1:2])

	if(nlayer > 1)
		bx = vcat(bx,  x[2:end-1], x[end-1]*ones(nz-3), 
		  x[end-2:-1:2], x[2]*ones(nz-4))

		bz = vcat(bz,  z[2]*ones(nx-2), z[3:end-1],
		  z[end-1]*ones(nx-3), z[end-2:-1:3])
	end
	if(nlayer > 2)
		bx = vcat(bx,  x[3:end-2], x[end-2]*ones(nz-5), 
		  x[end-3:-1:3], x[3]*ones(nz-6)
		  )
		bz = vcat(bz,  z[3]*ones(nx-4), z[4:end-2],
		  z[end-2]*ones(nx-5), z[end-3:-1:4]
		  )
	end
	isequal(length(bz), length(bx)) ? nothing : error("unequal dimensions")
	if(onlycount)
		return length(bz)
	else
		return bz, bx, length(bz)
	end
end


"""
Data type to represent 1D grid.
# Fields
* `x` : values
* `nx` : number of samples
* `δx` : sampling interval
"""
type M1D
	x::Array{Float64}
	nx::Int64
	δx::Float64
	"adding conditions that are to be false while construction"
	M1D(x, nx, δx) = 
		any([δx < 0.0, length(x) != nx]) ? 
			error("M1D construct") : new(x, nx, δx)
end

"Logical operation for `M1D`"
function M1D_isequal(grid1::M1D, grid2::M1D)
	return all([(isequal(getfield(grid1, name),getfield(grid2, name))) for name in fieldnames(M1D)])
end

"""
Construct 1-D grid based on nx or δx
"""
function M1D(xbeg::Float64, xend::Float64, nx::Int64)
	x = Array(linspace(xbeg, xend, nx))
	δx = length(x)==1 ? 0. : x[2]-x[1]
	return M1D(x, nx, δx)
end

function M1D(xbeg::Float64, xend::Float64, δx::Float64)
	x = [tt for tt in xbeg:δx:xend]
	δx = length(x)==1 ? 0. : x[2]-x[1]
	return M1D(x, size(x,1), δx)
end

"""
Grid with both positive and negative samples for a given lag.
Make sure that the number so samples is odd
"""
function M1D_lag(xlag::Float64, δx::Float64)
	x1=[tt for tt in 0.0:-δx:-xlag]
	x = vcat(flipdim(x1,1), -1.0.*x1[2:end])
	δx = length(x)==1 ? 0. : x[2]-x[1]
	isodd(length(x)) ? nothing : error("error in creating lag grid")
	return M1D(x, size(x,1), δx)
end

"""
1-D grid with a different sampling interval
* Not yet implemented for npow2 grids
"""
M1D_resamp(grid::M1D, δx::Float64) = M1D(grid.x[1], grid.x[end], δx)

"""
1-D grid which is has a different size
"""
function M1D_truncate(grid::M1D, xbeg::Float64, xend::Float64)
	ix1 = indmin((grid.x - xbeg).^2.0);
	ix2 = indmin((grid.x - xend).^2.0);
	ixmin = minimum([ix1, ix2]); 
	ixmax = maximum([ix1, ix2]); 
	x = grid.x[ixmin:ixmax];
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
