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


function interp_spray!(xin::Vector{Float64}, yin::AbstractVector{Float64},
					   xout::Vector{Float64}, yout::AbstractVector{Float64},
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
	ivec=zeros(Int64,np)

	if(attrib == :interp)
		# yout is only updated "within" bounds of yin 
		# get two indices and choose the inner most one
		ioutmin, ioutmax = indminn_inside(xout, xin[1], xin[end])

		yout[ioutmin:ioutmax] = zero(Float64);
		@simd for i in ioutmin:ioutmax
			indminn!(ivec,xin,xout[i]);
			interp_func(i, ivec, yout, xin,yin, xout[i])
		end
	elseif(attrib == :spray)
		# spary values of yin "within" bounds of yout
		# get 2 indices and choose the inner most one
		iinmin, iinmax = indminn_inside(xin,xout[1],xout[end])

		yout[:] = zero(Float64);
		@simd for i in iinmin:iinmax
			indminn!(ivec,xout,xin[i]);
			spray_func(ivec, yout, xout, xin[i], yin[i])
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
	ivec=zeros(Int64,np)

	if(attrib == :interp)
		# yout is only updated within bounds of yin 
		# get 2 indices and choose the inner most one
		i1outmin, i1outmax = indminn_inside(zout,zin[1],zin[end])
		i2outmin, i2outmax = indminn_inside(xout,xin[1],xin[end])
		yout[i1outmin:i1outmax, i2outmin:i2outmax] = zero(Float64);

		y_x=zeros(Float64, length(zin), length(i2outmin:i2outmax))
		# first along x
		@simd for ix in i2outmin:i2outmax
			@simd for iz in eachindex(zin)
				indminn!(ivec,xin, xout[ix]);
				iy=(ix-i2outmin)*length(zin)+iz
				yyin=view(yin, iz, :)
				interp_func(iy, ivec, y_x, xin, yyin, xout[ix])
			end
		end
		# then along z
		@simd for ix in i2outmin:i2outmax
			@simd for iz in i1outmin:i1outmax
				indminn!(ivec,zin,zout[iz]);
				iy=iz+(ix-1)*length(zout)
				yy_x=view(y_x, :, ix-i2outmin+1)
				interp_func(iy, ivec, yout, zin, yy_x, zout[iz])
			end
		end
	elseif(attrib == :spray)
		# spary values of yin only within bounds of yout
		# get 2 indices and choose the inner most one
		i1inmin, i1inmax = indminn_inside(zin,zout[1],zout[end])
		i2inmin, i2inmax = indminn_inside(xin,xout[1],xout[end])
		yout[:, :] = zero(Float64);
		y_z = zeros(Float64, length(zout), length(i2inmin:i2inmax))
		# first along z
		@simd for ix in i2inmin:i2inmax
			@simd for iz in i1inmin:i1inmax
				indminn!(ivec,zout,zin[iz]);
				yy_z=view(y_z,:,ix-i2inmin+1)
				spray_func(ivec, yy_z, zout, zin[iz], yin[iz, ix])
			end
		end
		# then along x
		@simd for ix in i2inmin:i2inmax
			@simd for iz in eachindex(zout)
				indminn!(ivec,xout,xin[ix]);
				yyout=view(yout,iz,:)
				spray_func(ivec, yyout, xout,  xin[ix], y_z[iz,ix-i2inmin+1])
			end
		end
	end
end # interp_spray!


"""
this subroutine interpolates or sprays bilinearly [one is adjoint of another]
interpolation returns y using Y[ivec[1]], Y[ivec[2]]
spraying returns Y[ivec[1]], Y[ivec[2]] using y
 
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
spraying returns Y[ivec[1]], y2, Y[ivec[3]], Y[ivec[4]] using y

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
		y[iy] =  (Float64(0.25) * (Y[ivec[1]] + Y[ivec[2]] + Y[ivec[3]] + Y[ivec[4]]))
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
