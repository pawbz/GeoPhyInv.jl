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
function indminn(x::AbstractVector{Float64}, val::Float64, n::Int64=1)
	# using enumerate to avoid indexing
	ivec = fill(0,n);
	xc = zeros(n);
	for inn in 1:n
		min_i = 0
		min_x = Inf
		for (i, xi) in enumerate(x)
			dist = abs(xi - val)
			if ((dist < min_x))
				min_x = dist
				min_i = i
				xc[]
			end
		end
		xc[inn] = x[min_i]
		x[min_i] = typemax(Float64)
		ivec[inn] = min_i
	end
	x[ivec] = xc
	return sort(ivec)
end

## slower version for large vectors
#function indminn(x::AbstractVector{Float64}, val::Float64, n::Int64)
#	xa =abs(x-val)
#	((n >= 1) & (n <= length(x))) ? nothing : error("invalid n")
#	ivec = [];
#	xc = [];
#	for i in 1:n
#		ii = indmin(xa);
#		push!(ivec, ii)
#		push!(xc, x[ii])
#		xa[ii] = typemax(Float64);
#	end
#	xa[ivec] = xc
#	return sort(ivec)
#end

function interp_spray!(xin::Vector{Float64}, yin::Vector{Float64},
					   xout::Vector{Float64}, yout::Vector{Float64},
					   attrib::Symbol,
					   Battrib::Symbol=:B1
					  )

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

	if(attrib == :interp)
		# yout is only updated within bounds of yin 
		ioutmin = indminn(xout,xin[1],1)[1]
		ioutmax = indminn(xout,xin[end],1)[1]

		yout[ioutmin:ioutmax] = zero(Float64);
		@simd for i in ioutmin:ioutmax
			ivec=indminn(xin,xout[i], np);
			interp_func(i, yout, xin[ivec]...,yin[ivec]..., xout[i])
		end
	elseif(attrib == :spray)
		# spary values of yin within bounds of yout
		iinmin = indminn(xin,xout[1],1)[1]
		iinmax = indminn(xin,xout[end],1)[1]

		yout[:] = zero(Float64);
		@simd for i in iinmin:iinmax
			ivec=indminn(xout,xin[i], np);
			spray_func(ivec..., yout, xout[ivec]..., xin[i], yin[i])
		end
	end
end # interp_spray!

function interp_spray!(xin::Array{Float64,1}, zin::Array{Float64,1}, yin::Array{Float64,2},
					   xout::Array{Float64,1}, zout::Array{Float64,1}, yout::Array{Float64,2},
					   attrib::Symbol, 
					   Battrib::Symbol=:B1
					  )
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

	if(attrib == :interp)
		# yout is only updated within bounds of yin 
		i1outmin = indminn(zout,zin[1],1)[1]
		i1outmax = indminn(zout,zin[end],1)[1]
		i2outmin = indminn(xout,xin[1],1)[1]
		i2outmax = indminn(xout,xin[end],1)[1]
		yout[i1outmin:i1outmax, i2outmin:i2outmax] = zero(Float64);

		y_x=zeros(Float64, length(zin), length(i2outmin:i2outmax))
		# first along x
		for ix in i2outmin:i2outmax
			@simd for iz in eachindex(zin)
				ivec=indminn(xin, xout[ix], np);
				iy=(ix-i2outmin)*length(zin)+iz
				interp_func(iy, y_x, xin[ivec]..., yin[iz, ivec]..., xout[ix])
			end
		end
		# then along z
		for ix in i2outmin:i2outmax
			@simd for iz in i1outmin:i1outmax
				ivec=indminn(zin,zout[iz], np);
				iy=iz+(ix-1)*length(zout)
				interp_func(iy, yout, zin[ivec]..., y_x[ivec, ix-i2outmin+1]..., zout[iz]...)
			end
		end
	elseif(attrib == :spray)
		# spary values of yin only within bounds of yout
		i1inmin = indminn(zin,zout[1],1)[1]
		i1inmax = indminn(zin,zout[end],1)[1]
		i2inmin = indminn(xin,xout[1],1)[1]
		i2inmax = indminn(xin,xout[end],1)[1]
		yout[:, :] = zero(Float64);
		y_z = zeros(Float64, length(zout), length(i2inmin:i2inmax))
		# first along z
		for ix in i2inmin:i2inmax
			@simd for iz in i1inmin:i1inmax
				ivec=indminn(zout,zin[iz], np);
				iy=ivec+(ix-i2inmin)*length(zout)
				spray_func(iy..., y_z, zout[ivec]..., zin[iz], yin[iz, ix])
			end
		end
		# then along x
		for ix in i2inmin:i2inmax
			@simd for iz in eachindex(zout)
				ivec=indminn(xout,xin[ix], np);
				iy=(ivec-1).*length(zout) + iz
				spray_func(iy..., yout, xout[ivec]...,  xin[ix], y_z[iz,ix-i2inmin+1])
			end
		end
	end
end # interp_spray!


"""
this subroutine interpolates or sprays bilinearly [one is adjoint of another]
interpolation returns y using y1, y2
spraying returns y1, y2 using y
 
                        +                      
                        |                      
    y1= f(x1)           |      y2= f(x2)   
      +-----------------x--------+             
                 y=f(x) |                      
                        +                      

bilinear interpolation
Reference: http://www.ajdesigner.com/phpinterpolation/bilinear_interpolation_equation.php
"""
function interp_B1_1D(iy, yV, x1, x2, y1, y2, x)
	denom = (x2 - x1)
	yV[iy] = ((y1*(x2-x) + y2*(x-x1)) / denom)
end # interpolate_spray_B1_1D

function spray_B1_1D( iy1, iy2, yV, x1, x2, x, y)
	denom = (x2 - x1)
	yV[iy1] += (y * (x2-x)/denom)
	yV[iy2] += (y * (x-x1)/denom)
end # interpolate_spray_B1_1D


"""
this subroutine interpolates or sprays 
using cubic bspline
interpolation returns y using y1, y2, y3, y4
spraying returns y1, y2, y3, y4 using y

                       +                      
                       |                      
   y1      y2          |    y3       y4
     +-------+---------x----+--------+             
                y=f(x) |                      
                       +                      
"""
function interp_B2_1D!(iy, yV, x1, x2, x3, x4, y1, y2, y3, y4, x)
	
	"some constants"
	h = (x2 - x1);
	hcube = h*h*h;
	hcubeI = hcube^(-1.0)

	if(h == zero(Float64)) then
		y[iy] =  (Float64(0.25) * (y1 + y2 + y3 + y4))
	else
		
		if((((x-x1)*(x-x2)) < zero(Float64)) | (x==x1)) 
			"left edge interpolate"
			c1 = (1.0 / 6.0 * (2.0 * h - (x - (x1-h)))^3) * hcubeI
			c4 = (1.0 / 6.0 * (2.0 * h - ((x4-h) - x))^3) * hcubeI
			c2 = ((2.0 / 3.0 * hcube - (0.5 * ((x2-h) - x)^2.0 * (2.0 * h - (x - (x2-h))) ) )) * hcubeI
			c3 = ((2.0 / 3.0 * hcube - (0.5 * ((x3-h) - x)^2.0 * (2.0 * h + (x - (x3-h))) ) )) * hcubeI
			yV[iy] = (y1 * (c2 + 2.0*c1) + y2 * (c3 - c1) + y3 * c4)
		elseif((((x-x3)*(x-x4)) < zero(Float64)) | (x==x4)) 
			"right edge interpolate"
			c1 = (1.0 / 6.0 * (2.0 * h - (x - (x1+h)))^3) * hcubeI
			c4 = (1.0 / 6.0 * (2.0 * h - ((x4+h) - x))^3) * hcubeI
			c2 = ((2.0 / 3.0 * hcube - (0.5 * ((x2+h) - x)^2.0 * (2.0 * h - (x - (x2+h))) ) )) * hcubeI
			c3 = ((2.0 / 3.0 * hcube - (0.5 * ((x3+h) - x)^2.0 * (2.0 * h + (x - (x3+h))) ) )) * hcubeI
			yV[iy] = (y2 * c1 + y3 * (c2 - c4) + y4 * (c3 + 2.0*c4))
		else
			c1 = (1.0 / 6.0 * (2.0 * h - (x - x1))^3) * hcubeI
			c4 = (1.0 / 6.0 * (2.0 * h - (x4 - x))^3) * hcubeI
			c2 = ((2.0 / 3.0 * hcube - (0.5 * (x2 - x)^2.0 * (2.0 * h - (x - x2)) ) )) * hcubeI
			c3 = ((2.0 / 3.0 * hcube - (0.5 * (x3 - x)^2.0 * (2.0 * h + (x - x3)) ) )) * hcubeI
			"center interpolate"
			yV[iy] = ((y1 * c1 + y2 * c2 + y3 * c3 + y4 * c4))
		end
	end
	return yV[iy]

end # interp_B2_1D

function spray_B2_1D!(iy1, iy2, iy3, iy4, yV, x1, x2, x3, x4, x, y)

	"some constants"
	h = (x2 - x1);
	hcube = h*h*h;
	hcubeI = hcube^(-1.0)

	if(h == zero(Float64)) then
		yV[iy1] += (y)
		yV[iy2] += (y) 
		yV[iy3] += (y) 
		yV[iy4] += (y) 
	else
		if((((x-x1)*(x-x2)) < zero(Float64)) | (x==x1)) 
			"left edge spray"
			c1 = (1.0 / 6.0 * (2.0 * h - (x - (x1-h)))^3) * hcubeI
			c4 = (1.0 / 6.0 * (2.0 * h - ((x4-h) - x))^3) * hcubeI
			c2 = ((2.0 / 3.0 * hcube - (0.5 * ((x2-h) - x)^2.0 * (2.0 * h - (x - (x2-h))) ) )) * hcubeI
			c3 = ((2.0 / 3.0 * hcube - (0.5 * ((x3-h) - x)^2.0 * (2.0 * h + (x - (x3-h))) ) )) * hcubeI

			yV[iy1] += (y * (c2 + 2.0*c1))
			yV[iy2] += (y * (c3 - c1))
			yV[iy3] += (y * c4)
		elseif((((x-x3)*(x-x4)) < zero(Float64)) | (x==x4)) 
			"right edge spray"        
			c1 = (1.0 / 6.0 * (2.0 * h - (x - (x1+h)))^3) * hcubeI
			c4 = (1.0 / 6.0 * (2.0 * h - ((x4+h) - x))^3) * hcubeI
			c2 = ((2.0 / 3.0 * hcube - (0.5 * ((x2+h) - x)^2.0 * (2.0 * h - (x - (x2+h))) ) )) * hcubeI
			c3 = ((2.0 / 3.0 * hcube - (0.5 * ((x3+h) - x)^2.0 * (2.0 * h + (x - (x3+h))) ) )) * hcubeI

			yV[iy2] += y * c1 
			yV[iy3] += y * (c2 - c4) 
			yV[iy4] += y * (c3+ 2.0 * c4) 
		else
			"center spray"
			c1 = (1.0 / 6.0 * (2.0 * h - (x - x1))^3) * hcubeI
			c4 = (1.0 / 6.0 * (2.0 * h - (x4 - x))^3) * hcubeI
			c2 = ((2.0 / 3.0 * hcube - (0.5 * (x2 - x)^2.0 * (2.0 * h - (x - x2)) ) )) * hcubeI
			c3 = ((2.0 / 3.0 * hcube - (0.5 * (x3 - x)^2.0 * (2.0 * h + (x - x3)) ) )) * hcubeI
			yV[iy1] += (y * c1)
			yV[iy2] += (y * c2)
			yV[iy3] += (y * c3) 
			yV[iy4] += (y * c4) 
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
