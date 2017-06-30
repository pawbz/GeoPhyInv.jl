__precompile__()

module Interpolation
"""
## TODO:
* add dimension checks to interp_spray!
Reference: https://www.ibiblio.org/e-notes/Splines/bezier.html

Return n indices in order
Cannot find a julia method which does, this.
If a faster method is found, replace it later.
"""
function indminn{T<:Real}(x::Vector{T}, n::Int64)
	xc = copy(x);
	((n >= 1) & (n <= length(x))) ? nothing : error("invalid n")
	ivec = [];
	for i in 1:n
		ii = indmin(xc);
		xc[ii] = typemax(T);
		push!(ivec, ii)
	end
	return sort(ivec)
end

function interp_spray!{T1<:Real, T2<:Real}(xin::Vector{T1}, yin::Vector{T2},
					   xout::Vector{T1}, yout::Vector{T2},
					   attrib::Symbol,
					   Battrib::Symbol=:B1
					  )

	if(Battrib == :B1)
		np=2;
		interp_func = interp_B1_1D
		spray_func = spray_B1_1D
	elseif(Battrib == :B2)
		np=4;
		interp_func = interp_B2_1D
		spray_func = spray_B2_1D
	else
		error("invalid attrib")
	end

	yout[:] = zero(T2);
	if(attrib == :interp)
		for i in eachindex(xout)
			ivec=indminn(abs(xin-xout[i]), np);
			yout[i] = interp_func(xin[ivec], yin[ivec], xout[i])
		end
	elseif(attrib == :spray)
		for i in eachindex(xin)
			ivec=indminn(abs(xout-xin[i]), np);
			yout[ivec] += spray_func(xout[ivec], xin[i], yin[i])
		end
	end
end # interp_spray!

function interp_spray!{T1<:Real, T2<:Real}(xin::Array{T1,1}, zin::Array{T1,1}, yin::Array{T2,2},
					   xout::Array{T1,1}, zout::Array{T1,1}, yout::Array{T2,2},
					   attrib::Symbol, 
					   Battrib::Symbol=:B1
					  )
	if(Battrib == :B1)
		np=2;
		interp_func = interp_B1_1D
		spray_func = spray_B1_1D
	elseif(Battrib == :B2)
		np=4;
		interp_func = interp_B2_1D
		spray_func = spray_B2_1D
	else
		error("invalid attrib")
	end

	yout[:,:] = zero(T2);
	if(attrib == :interp)
		y_x=zeros(T2, length(zin), length(xout))
		# first along x
		for iz in eachindex(zin)
			for ix in eachindex(xout)
				ivec=indminn(abs(xin-xout[ix]), np);
				y_x[iz,ix] = interp_func(xin[ivec], yin[iz,ivec], xout[ix])
			end
		end
		# then along z
		for ix in eachindex(xout)
			for iz in eachindex(zout)
				ivec=indminn(abs(zin-zout[iz]), np);
				yout[iz,ix] = interp_func(zin[ivec], y_x[ivec,ix], zout[iz])
			end
		end
	elseif(attrib == :spray)
		y_z = zeros(T2, length(zout), length(xin))
		# first along z
		for ix in eachindex(xin)
			for iz in eachindex(zin)
				ivec=indminn(abs(zout-zin[iz]), np);
				y_z[ivec,ix] += spray_func(zout[ivec], zin[iz], yin[iz,ix])
			end
		end
		# then along x
		for iz in eachindex(zout)
			for ix in eachindex(xin)
				ivec=indminn(abs(xout-xin[ix]), np);
				yout[iz,ivec] += spray_func(xout[ivec], xin[ix], y_z[iz,ix])
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
function interp_B1_1D{T1<:Real,T2<:Real}(xV::Vector{T1}, yV::Vector{T2}, x::T1)
	denom = (xV[2] - xV[1])
	return ((yV[1]*(xV[2]-x) + yV[2]*(x-xV[1])) / denom)
end # interpolate_spray_B1_1D

function spray_B1_1D{T1<:Real,T2<:Real}(xV::Vector{T1}, x::T1, y::T2)
	denom = (xV[2] - xV[1])
	yV = zeros(T2,length(xV))
	if((xV[2]-x) == zero(T1)) 
	else
		yV[1] = (y * (xV[2]-x)/denom)
	end
	if((x-xV[1]) == zero(T1))
		yV[2] = (zero(T2));
	else
		yV[2] = (y * (x-xV[1])/denom)
	end
	return yV
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
function interp_B2_1D{T1<:Real,T2<:Real}(xV::Vector{T1}, yV::Vector{T2}, x::T1)
	
	"some constants"
	h = (xV[2] - xV[1]);
	hcube = h*h*h;
	hcubeI = hcube^(-1.0)

	if(h == zero(T1)) then
		return (T(0.25) * (yV[1] + yV[2] + yV[3] + yV[4]))
	else
		
		if((((x-xV[1])*(x-xV[2])) < zero(T1)) | (x==xV[1])) 
			"left edge interpolate"
			c1 = (1.0 / 6.0 * (2.0 * h - (x - (xV[1]-h)))^3) * hcubeI
			c4 = (1.0 / 6.0 * (2.0 * h - ((xV[4]-h) - x))^3) * hcubeI
			c2 = ((2.0 / 3.0 * hcube - (0.5 * ((xV[2]-h) - x)^2.0 * (2.0 * h - (x - (xV[2]-h))) ) )) * hcubeI
			c3 = ((2.0 / 3.0 * hcube - (0.5 * ((xV[3]-h) - x)^2.0 * (2.0 * h + (x - (xV[3]-h))) ) )) * hcubeI
			return (yV[1] * (c2 + 2.0*c1) + yV[2] * (c3 - c1) + yV[3] * c4)
		elseif((((x-xV[3])*(x-xV[4])) < zero(T1)) | (x==xV[4])) 
			"right edge interpolate"
			c1 = (1.0 / 6.0 * (2.0 * h - (x - (xV[1]+h)))^3) * hcubeI
			c4 = (1.0 / 6.0 * (2.0 * h - ((xV[4]+h) - x))^3) * hcubeI
			c2 = ((2.0 / 3.0 * hcube - (0.5 * ((xV[2]+h) - x)^2.0 * (2.0 * h - (x - (xV[2]+h))) ) )) * hcubeI
			c3 = ((2.0 / 3.0 * hcube - (0.5 * ((xV[3]+h) - x)^2.0 * (2.0 * h + (x - (xV[3]+h))) ) )) * hcubeI
			return (yV[2] * c1 + yV[3] * (c2 - c4) + yV[4] * (c3 + 2.0*c4))
		else
			c1 = (1.0 / 6.0 * (2.0 * h - (x - xV[1]))^3) * hcubeI
			c4 = (1.0 / 6.0 * (2.0 * h - (xV[4] - x))^3) * hcubeI
			c2 = ((2.0 / 3.0 * hcube - (0.5 * (xV[2] - x)^2.0 * (2.0 * h - (x - xV[2])) ) )) * hcubeI
			c3 = ((2.0 / 3.0 * hcube - (0.5 * (xV[3] - x)^2.0 * (2.0 * h + (x - xV[3])) ) )) * hcubeI
			"center interpolate"
			return ((yV[1] * c1 + yV[2] * c2 + yV[3] * c3 + yV[4] * c4))
		end
	end

end # interp_B2_1D

function spray_B2_1D{T1<:Real,T2<:Real}(xV::Vector{T1}, x::T1, y::T2)

	"some constants"
	h = (xV[2] - xV[1]);
	hcube = h*h*h;
	hcubeI = hcube^(-1.0)

	yV = zeros(T2,length(xV))
	if(h == zero(T1)) then
		yV[1] = copy(y)
		yV[2] = copy(y) 
		yV[3] = copy(y) 
		yV[4] = copy(y) 
		return yV
	else
		if((((x-xV[1])*(x-xV[2])) < zero(T1)) | (x==xV[1])) 
			"left edge spray"
			c1 = (1.0 / 6.0 * (2.0 * h - (x - (xV[1]-h)))^3) * hcubeI
			c4 = (1.0 / 6.0 * (2.0 * h - ((xV[4]-h) - x))^3) * hcubeI
			c2 = ((2.0 / 3.0 * hcube - (0.5 * ((xV[2]-h) - x)^2.0 * (2.0 * h - (x - (xV[2]-h))) ) )) * hcubeI
			c3 = ((2.0 / 3.0 * hcube - (0.5 * ((xV[3]-h) - x)^2.0 * (2.0 * h + (x - (xV[3]-h))) ) )) * hcubeI

			yV[1] = copy(y * (c2 + 2.0*c1))
			yV[2] = copy(y * (c3 - c1))
			yV[3] = copy(y * c4)
			return yV
		elseif((((x-xV[3])*(x-xV[4])) < zero(T1)) | (x==xV[4])) 
			"right edge spray"        
			c1 = (1.0 / 6.0 * (2.0 * h - (x - (xV[1]+h)))^3) * hcubeI
			c4 = (1.0 / 6.0 * (2.0 * h - ((xV[4]+h) - x))^3) * hcubeI
			c2 = ((2.0 / 3.0 * hcube - (0.5 * ((xV[2]+h) - x)^2.0 * (2.0 * h - (x - (xV[2]+h))) ) )) * hcubeI
			c3 = ((2.0 / 3.0 * hcube - (0.5 * ((xV[3]+h) - x)^2.0 * (2.0 * h + (x - (xV[3]+h))) ) )) * hcubeI

			yV[2] = y * c1 
			yV[3] = y * (c2 - c4) 
			yV[4] = y * (c3+ 2.0 * c4) 
			return yV
		else
			"center spray"
     			c1 = (1.0 / 6.0 * (2.0 * h - (x - xV[1]))^3) * hcubeI
			c4 = (1.0 / 6.0 * (2.0 * h - (xV[4] - x))^3) * hcubeI
			c2 = ((2.0 / 3.0 * hcube - (0.5 * (xV[2] - x)^2.0 * (2.0 * h - (x - xV[2])) ) )) * hcubeI
			c3 = ((2.0 / 3.0 * hcube - (0.5 * (xV[3] - x)^2.0 * (2.0 * h + (x - xV[3])) ) )) * hcubeI
			yV[1] = copy(y * c1)
			yV[2] = copy(y * c2)
			yV[3] = copy(y * c3) 
			yV[4] = copy(y * c4) 
			return yV
		end
	end
end # interpolate_spray_B3_1D


"""
get source spray weights
"""
function get_spray_weights!(weights::AbstractArray{Float64}, denomI::AbstractArray{Float64}, 
			   ix1::AbstractArray{Int64}, ix2::AbstractArray{Int64}, iz1::AbstractArray{Int64}, iz2::AbstractArray{Int64}, 
			   mesh_x::Vector{Float64}, mesh_z::Vector{Float64}, xval::Float64, zval::Float64)


	ix1[1], ix2[1] = indminn(abs(mesh_x-xval),2)
	iz1[1], iz2[1] = indminn(abs(mesh_z-zval),2)
	
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

	ix1[1], ix2[1]=indminn(abs(mesh_x-xval),2)
	iz1[1], iz2[1]=indminn(abs(mesh_z-zval),2)

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
