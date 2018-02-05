__precompile__()

"""
## TODO:
* add dimension checks to interp_spray!
Reference: https://www.ibiblio.org/e-notes/Splines/bezier.html
"""
module Interpolation

"""
Return n indices in order
Cannot find a julia method which does, this.
If a faster method is found, replace it later.
"""
function indminn!(ivec, x::AbstractVector{Float64}, val::Float64)
	# using enumerate to avoid indexing
	n=length(ivec)
	ivec[:]=zero(eltype(ivec))
	for inn in 1:n
		ivec[inn] = indminimum(x, val, ivec)
	end
	sort!(ivec)
	return ivec
end

"""
minimum index using imask
"""
function indminimum(x, val, imask=[])
	min_i=0
	min_x=typemax(Float64)
	for (i, xi) in enumerate(x)
		dist = abs(xi - val)
		if ((dist < min_x) && (i ∉ imask))
			min_x = dist
			min_i = i
		end
	end
	return min_i
end

function indminn(x::AbstractVector{Float64}, val::Float64, n::Int64=1)
	ivec=fill(0,n)
	indminn!(ivec,x,val)
	return ivec
end

"return index such that "
function indminn_inside(x, valbeg, valend)
	iibeg1 = indminimum(x, valbeg)
	iibeg2 = indminimum(x, valbeg,iibeg1)
	iiend1 = indminimum(x, valend)
	iiend2 = indminimum(x, valend,iiend1)

	if((valbeg-x[iibeg1])*(valend-x[iibeg1]) <= 0.0)
		iibeg=iibeg1
	else
		iibeg=iibeg2
	end
	if((valbeg-x[iiend1])*(valend-x[iiend1])<=0.0)
		iiend=iiend1
	else
		iiend=iiend2
	end

	return iibeg, iiend
end

"""
N is either one or two
* `x` coordinates before interpolation
* `xi` coordinates after interpolation
"""
type Param{T<:Real}
	x::Array{Array{T,1},1}
	xi::Array{Array{T,1},1}
	np::Int
	nd::Int
	nx::Vector{Int}
	nxi::Vector{Int}
	ivecinterp::Vector{Matrix{Int64}}
	ivecspray::Vector{Matrix{Int64}}
	iximin::Vector{Int64}
	iximax::Vector{Int64}
	ixmin::Vector{Int64}
	ixmax::Vector{Int64}
	interp_func::Function
	spray_func::Function
	Battrib::Symbol # select spline 
	y_x::Matrix{T}
	y_z::Matrix{T}
end


function Param(x::Array{Vector{Float64}}, xi::Array{Vector{Float64}}, Battrib::Symbol=:B1)

	if(Battrib == :B1)
		np=2;
		interp_func = interp_B1_1D 
		spray_func = spray_B1_1D 
	elseif(Battrib == :B2)
		np=4;
		interp_func = interp_B2_1D!
		spray_func = spray_B2_1D!
	else
		error("invalid attrib")
	end
	!(length(x)==length(xi)) && error("incorrect dimensions x")
	nd=length(x)
	nx=[length(x[id]) for id in 1:nd]
	nxi=[length(xi[id]) for id in 1:nd]

	ivecinterp=[zeros(Int64,np,nxi[id]) for id in 1:nd]
	ivecspray=[zeros(Int64,np,nx[id]) for id in 1:nd]

	iximin=zeros(Int64,nd)
	iximax=zeros(Int64,nd)
	ixmin=zeros(Int64,nd)
	ixmax=zeros(Int64,nd)

	for id in 1:nd
		# yi is only updated "within" bounds of y 
		# get two indices and choose the inner most one
		iximin[id], iximax[id] = indminn_inside(xi[id], x[id][1], x[id][end])
		# spary values of y "within" bounds of yi
		# get 2 indices and choose the inner most one
		ixmin[id], ixmax[id] = indminn_inside(x[id],xi[id][1],xi[id][end])

		# index interp
		for i in 1:nxi[id]
			ivec=ivecinterp[id]
			ivecc=view(ivec, :, i)
			xn=x[id]
			xit=xi[id]
			indminn!(ivecc,xn,xit[i])
		end
		# index spray
		for i in 1:nx[id]
			ivec=ivecspray[id]
			ivecc=view(ivec, :, i)
			xn=x[id]
			xit=xi[id]
			indminn!(ivecc,xit,xn[i])
		end

	end

	# some other allocations for 2D
	if(nd==2)
		y_x=zeros(Float64, length(x[2]), length(iximin[1]:iximax[1]))
		y_z = zeros(Float64, length(x[2]), length(iximin[1]:iximax[1]))
	else
		y_x=zeros(1,1)
		y_z=zeros(1,1)
	end
	return Param(
	x,	xi,	np,	nd,	nx,	nxi,	ivecinterp,	ivecspray,	iximin,	iximax,	ixmin,
	ixmax,	interp_func,	spray_func,	Battrib, y_x, y_z)


end # Param


function interp_spray!(y::AbstractVector{Float64}, yi::AbstractVector{Float64} , pa::Param{Float64}, attrib)

	if(attrib == :interp)
		# yi is only updated "within" bounds of y 
		for i in pa.iximin[1] : pa.iximax[1]
			yi[i] = zero(Float64)
		end
		ivec=pa.ivecinterp[1]
		@simd for i in pa.iximin[1]:pa.iximax[1]
			ivecc=view(ivec,:,i)
			pa.interp_func(i, ivecc, yi, pa.x[1], y, pa.xi[1][i])
		end
	elseif(attrib == :spray)
		for i in eachindex(y)
			y[i] = zero(Float64)
		end
		# spary values of yi "within" bounds of y
		ivec=pa.ivecinterp[1]
		@simd for i in pa.iximin[1]:pa.iximax[1]
			ivecc=view(ivec,:,i)
			pa.spray_func(ivecc, y, pa.x[1], pa.xi[1][i], yi[i])
		end
	end
end # interp_spray!

"""
This allocates memory, need to be fixed
"""
function interp_spray!(y::Matrix{Float64}, yi::Matrix{Float64} , pa::Param{Float64}, attrib)

	if(attrib == :interp)
		for i1 in pa.iximin[2]:pa.iximax[2], i2 in pa.iximin[1]:pa.iximax[1]
			yi[i1, i2] = zero(Float64);
		end
		pa.y_x[:]=zero(Float64)

		# first along x
		ivec=pa.ivecinterp[1]
		@simd for ix in pa.iximin[1]:pa.iximax[1]
			@simd for iz in eachindex(pa.x[2])
				ivecc=view(ivec,:,ix)
				iy=(ix-pa.iximin[1])*length(pa.x[2])+iz
				yy=view(y, iz, :)
				pa.interp_func(iy, ivecc, pa.y_x, pa.x[1], yy, pa.xi[1][ix])
			end
		end
		# then along z
		ivec=pa.ivecinterp[2]
		@simd for ix in pa.iximin[1]:pa.iximax[1]
			@simd for iz in pa.iximin[2]:pa.iximax[2]
				ivecc=view(ivec,:,iz)
				iy=iz+(ix-1)*length(pa.xi[2])
				yy_x=view(pa.y_x, :, ix-pa.iximin[1]+1)
				pa.interp_func(iy, ivecc, yi, pa.x[2], yy_x, pa.xi[2][iz])
			end
		end
	elseif(attrib == :spray)
		# spary values of y only within bounds of yi
		# get 2 indices and choose the inner most one
		for i in eachindex(y)
			y[i]=zero(Float64)
		end
		pa.y_z[:]=zero(Float64)
		# first along z
		ivec=pa.ivecinterp[2]
		@simd for ix in pa.iximin[1]:pa.iximax[1]
			@simd for iz in pa.iximin[2]:pa.iximax[2]
				ivecc=view(ivec,:,iz)
				yy_z=view(pa.y_z,:,ix-pa.iximin[1]+1)
				pa.spray_func(ivecc, yy_z, pa.x[2], pa.xi[2][iz], yi[iz, ix])
			end
		end
		# then along x
		ivec=pa.ivecinterp[1]
		@simd for ix in  pa.iximin[1]:pa.iximax[1]
			@simd for iz in eachindex(pa.x[2])
				ivecc=view(ivec,:,ix)
				yy=view(y,iz,:)
				pa.spray_func(ivecc, yy, pa.x[1],  pa.xi[1][ix], pa.y_z[iz,ix-pa.iximin[1]+1])
			end
		end
	end
end # interp_spray!


"""
this subroutine interpolates or sprays bilinearly [one is adjoint of another]
interpolation returns y using Y[ivec[1]], Y[ivec[2]]
sprayg returns Y[ivec[1]], Y[ivec[2]] using y
 
                        +                      
                        |                      
    Y[ivec[1]]= f(X[ivec[1]])           |      Y[ivec[2]]= f(X[ivec[2]])   
      +-----------------x--------+             
                 y=f(x) |                      
                        +                      

bilinear interpolation
Reference: http://www.ajdesigner.com/phpinterpolation/bilinear_interpolation_equation.php
"""
function interp_B1_1D(iy, ivec, yV, X, Y, x)
	@inbounds yV[iy] = ((Y[ivec[1]]*(X[ivec[2]]-x) + Y[ivec[2]]*(x-X[ivec[1]])) * inv(X[ivec[2]]-X[ivec[1]]))
end # interpolate_spray_B1_1D

function spray_B1_1D(ivec, yV, X, x, y)
	denom = (X[ivec[2]] - X[ivec[1]])
	@inbounds yV[ivec[1]] += (y * (X[ivec[2]]-x)/denom)
	@inbounds yV[ivec[2]] += (y * (x-X[ivec[1]])/denom)
end # interpolate_spray_B1_1D


"""
this subroutine interpolates or sprays 
using cubic bspline
interpolation returns y using Y[ivec[1]], y2, Y[ivec[3]], Y[ivec[4]]
sprayg returns Y[ivec[1]], y2, Y[ivec[3]], Y[ivec[4]] using y

                       +                      
                       |                      
   Y[ivec[1]]      y2          |    Y[ivec[3]]       Y[ivec[4]]
     +-------+---------x----+--------+             
                y=f(x) |                      
                       +                      
"""
function interp_B2_1D!(iy, ivec, yV, X,Y,x)
	
	"some constants"
	h = (X[ivec[2]] - X[ivec[1]]);
	hcube = h*h*h;
	hcubeI = hcube^(-1.0)

	if(h == zero(Float64)) then
		yV[iy] =  (Float64(0.25) * (Y[ivec[1]] + Y[ivec[2]] + Y[ivec[3]] + Y[ivec[4]]))
	else
		
		if((((x-X[ivec[1]])*(x-X[ivec[2]])) < zero(Float64)) | (x==X[ivec[1]])) 
			"left edge interpolate"
			c1 = (1.0 / 6.0 * (2.0 * h - (x - (X[ivec[1]]-h)))^3) * hcubeI
			c4 = (1.0 / 6.0 * (2.0 * h - ((X[ivec[4]]-h) - x))^3) * hcubeI
			c2 = ((2.0 / 3.0 * hcube - (0.5 * ((X[ivec[2]]-h) - x)^2.0 * (2.0 * h - (x - (X[ivec[2]]-h))) ) )) * hcubeI
			c3 = ((2.0 / 3.0 * hcube - (0.5 * ((X[ivec[3]]-h) - x)^2.0 * (2.0 * h + (x - (X[ivec[3]]-h))) ) )) * hcubeI
			@inbounds yV[iy] = (Y[ivec[1]] * (c2 + 2.0*c1) + Y[ivec[2]] * (c3 - c1) + Y[ivec[3]] * c4)
		elseif((((x-X[ivec[3]])*(x-X[ivec[4]])) < zero(Float64)) | (x==X[ivec[4]])) 
			"right edge interpolate"
			c1 = (1.0 / 6.0 * (2.0 * h - (x - (X[ivec[1]]+h)))^3) * hcubeI
			c4 = (1.0 / 6.0 * (2.0 * h - ((X[ivec[4]]+h) - x))^3) * hcubeI
			c2 = ((2.0 / 3.0 * hcube - (0.5 * ((X[ivec[2]]+h) - x)^2.0 * (2.0 * h - (x - (X[ivec[2]]+h))) ) )) * hcubeI
			c3 = ((2.0 / 3.0 * hcube - (0.5 * ((X[ivec[3]]+h) - x)^2.0 * (2.0 * h + (x - (X[ivec[3]]+h))) ) )) * hcubeI
			@inbounds yV[iy] = (Y[ivec[2]] * c1 + Y[ivec[3]] * (c2 - c4) + Y[ivec[4]] * (c3 + 2.0*c4))
		else
			c1 = (1.0 / 6.0 * (2.0 * h - (x - X[ivec[1]]))^3) * hcubeI
			c4 = (1.0 / 6.0 * (2.0 * h - (X[ivec[4]] - x))^3) * hcubeI
			c2 = ((2.0 / 3.0 * hcube - (0.5 * (X[ivec[2]] - x)^2.0 * (2.0 * h - (x - X[ivec[2]])) ) )) * hcubeI
			c3 = ((2.0 / 3.0 * hcube - (0.5 * (X[ivec[3]] - x)^2.0 * (2.0 * h + (x - X[ivec[3]])) ) )) * hcubeI
			"center interpolate"
			@inbounds yV[iy] = ((Y[ivec[1]] * c1 + Y[ivec[2]] * c2 + Y[ivec[3]] * c3 + Y[ivec[4]] * c4))
		end
	end
	return yV[iy]

end # interp_B2_1D

function spray_B2_1D!(ivec, yV, X, x, y)

	"some constants"
	h = (X[ivec[2]] - X[ivec[1]]);
	hcube = h*h*h;
	hcubeI = hcube^(-1.0)

	if(h == zero(Float64)) then
		@inbounds yV[ivec[1]] += (y)
		@inbounds yV[ivec[2]] += (y) 
		@inbounds yV[ivec[3]] += (y) 
		@inbounds yV[ivec[4]] += (y) 
	else
		if((((x-X[ivec[1]])*(x-X[ivec[2]])) < zero(Float64)) | (x==X[ivec[1]])) 
			"left edge spray"
			c1 = (1.0 / 6.0 * (2.0 * h - (x - (X[ivec[1]]-h)))^3) * hcubeI
			c4 = (1.0 / 6.0 * (2.0 * h - ((X[ivec[4]]-h) - x))^3) * hcubeI
			c2 = ((2.0 / 3.0 * hcube - (0.5 * ((X[ivec[2]]-h) - x)^2.0 * (2.0 * h - (x - (X[ivec[2]]-h))) ) )) * hcubeI
			c3 = ((2.0 / 3.0 * hcube - (0.5 * ((X[ivec[3]]-h) - x)^2.0 * (2.0 * h + (x - (X[ivec[3]]-h))) ) )) * hcubeI

			@inbounds yV[ivec[1]] += (y * (c2 + 2.0*c1))
			@inbounds yV[ivec[2]] += (y * (c3 - c1))
			@inbounds yV[ivec[3]] += (y * c4)
		elseif((((x-X[ivec[3]])*(x-X[ivec[4]])) < zero(Float64)) | (x==X[ivec[4]])) 
			"right edge spray"        
			c1 = (1.0 / 6.0 * (2.0 * h - (x - (X[ivec[1]]+h)))^3) * hcubeI
			c4 = (1.0 / 6.0 * (2.0 * h - ((X[ivec[4]]+h) - x))^3) * hcubeI
			c2 = ((2.0 / 3.0 * hcube - (0.5 * ((X[ivec[2]]+h) - x)^2.0 * (2.0 * h - (x - (X[ivec[2]]+h))) ) )) * hcubeI
			c3 = ((2.0 / 3.0 * hcube - (0.5 * ((X[ivec[3]]+h) - x)^2.0 * (2.0 * h + (x - (X[ivec[3]]+h))) ) )) * hcubeI

			@inbounds yV[ivec[2]] += y * c1 
			@inbounds yV[ivec[3]] += y * (c2 - c4) 
			@inbounds yV[ivec[4]] += y * (c3+ 2.0 * c4) 
		else
			"center spray"
			c1 = (1.0 / 6.0 * (2.0 * h - (x - X[ivec[1]]))^3) * hcubeI
			c4 = (1.0 / 6.0 * (2.0 * h - (X[ivec[4]] - x))^3) * hcubeI
			c2 = ((2.0 / 3.0 * hcube - (0.5 * (X[ivec[2]] - x)^2.0 * (2.0 * h - (x - X[ivec[2]])) ) )) * hcubeI
			c3 = ((2.0 / 3.0 * hcube - (0.5 * (X[ivec[3]] - x)^2.0 * (2.0 * h + (x - X[ivec[3]])) ) )) * hcubeI
			@inbounds yV[ivec[1]] += (y * c1)
			@inbounds yV[ivec[2]] += (y * c2)
			@inbounds yV[ivec[3]] += (y * c3) 
			@inbounds yV[ivec[4]] += (y * c4) 
		end
	end
end # interpolate_spray_B3_1D


"""
get source spray weights
"""
function get_spray_weights!(weights::AbstractArray{Float64}, denomI::AbstractArray{Float64}, 
			   ix1::AbstractArray{Int64}, ix2::AbstractArray{Int64}, iz1::AbstractArray{Int64}, iz2::AbstractArray{Int64}, 
			   mesh_x::Vector{Float64}, mesh_z::Vector{Float64}, xval::Float64, zval::Float64)


	ix1[1], ix2[1] = indminn(mesh_x,xval,2)
	iz1[1], iz2[1] = indminn(mesh_z,zval,2)
	
	denomI[1] = ((mesh_x[ix2[1]] - mesh_x[ix1[1]])*(mesh_z[iz2[1]] - mesh_z[iz1[1]]))^(-1.e0)
	"for iz1, ix1"
	weights[1] = 	(mesh_x[ix2[1]]-	xval)*	(mesh_z[iz2[1]]-	zval)*	denomI[1] 
	"for iz1, ix2"
	weights[2] =   (xval-	mesh_x[ix1[1]])*	(mesh_z[iz2[1]]-	zval)*	denomI[1]
	"for iz2, ix1"
	weights[3] = 	(mesh_x[ix2[1]]-	xval)*	(zval-	mesh_z[iz1[1]])*	denomI[1]
	"for iz2, ix2" 
	weights[4] = 	(xval-	mesh_x[ix1[1]])*	(zval-	mesh_z[iz1[1]])*	denomI[1]
	return weights, denomI, ix1, ix2, iz1, iz2
end


"""
get interpolation weights
"""
function get_interpolate_weights!(weights::AbstractArray{Float64}, denomI::AbstractArray{Float64}, 
			   ix1::AbstractArray{Int64}, ix2::AbstractArray{Int64}, iz1::AbstractArray{Int64}, iz2::AbstractArray{Int64}, 
			   mesh_x::Vector{Float64}, mesh_z::Vector{Float64}, xval::Float64, zval::Float64)

	ix1[1], ix2[1] = indminn(mesh_x,xval,2)
	iz1[1], iz2[1] = indminn(mesh_z,zval,2)

	denomI[1] = ((mesh_x[ix2[1]] -mesh_x[ix1[1]])*(mesh_z[iz2[1]] - mesh_z[iz1[1]]))^(-1.e0)

	"iz1, ix1"
	weights[1] =	((mesh_x[ix2[1]]-	xval)*	(mesh_z[iz2[1]]-zval))*denomI[1]
	"iz1, ix2"
	weights[2] = ((xval-	mesh_x[ix1[1]])*(mesh_z[iz2[1]]-zval))*denomI[1]
	"recz2, ix1"
	weights[3] =	((mesh_x[ix2[1]]-	xval)*(zval-mesh_z[iz1[1]]))*denomI[1]
	"iz2, ix2"
	weights[4] =	((xval-	mesh_x[ix1[1]])*	(zval-	mesh_z[iz1[1]]))*	denomI[1]

end


end # module
