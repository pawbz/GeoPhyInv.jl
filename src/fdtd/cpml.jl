"""
Generate coefficients related to PML boundaries, e.g., damping profiles.
This function outputs a_x, b_x, and k_x variables in eq. 25-26 from Komatitsch 2007 (Geophysics).

	TODO
* At the moment, K=1, and is dummy, unlike in the case of EM when it was useful. Removing K and K_inv might save us memory when using GPUs.
* Restrict the sizes of these coefficients only to the PML region while still writing optimized loops?
"""
function pml_coeff(mgrid, # 1D grid
	        flags::Vector{Bool},
		δt::Float64, 
		np::Int64, 
		velmax::Float64, velmin::Float64, 
		freqpeak::Float64,
		)

	δx=step(mgrid); nx=length(mgrid)

	NPOWER = 2.e0
	K_MAX_PML = 1.e0 # from Gedney page 8.11
	ALPHA_MAX_PML = pi * freqpeak # from Festa and Vilotte

	# thickness of the PML layer in meters
	thickness_PML_x = np * δx

	"reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf"
	Rcoef = 0.001e0

	#! check that NPOWER is okay
	#if(NPOWER < 1) stop "NPOWER must be greater than 1"

	# compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
	d0_x = - (NPOWER + 1) * (velmax+velmin)/2.0 * log(Rcoef) / (2.e0 * thickness_PML_x)


	d_x=Data.Array(zeros(nx)); d_x_half=zero(d_x)
	k_x=zero(d_x);	k_xI=zero(d_x); fill!(k_x,1); fill!(k_xI,1)
	k_x_half=zero(d_x);	k_x_halfI=zero(d_x); fill!(k_x_half,1); fill!(k_x_halfI,1)
	alpha_x=zero(d_x);	alpha_x_half=zero(d_x)
	a_x=zero(d_x);	a_x_half=zero(d_x)
	b_x=zero(d_x);	b_x_half=zero(d_x)


	"damping in the X direction"
	"origin of the PML layer (position of right edge minus thickness, in meters)"
	xoriginleft = thickness_PML_x
	xoriginright = (nx-1)*δx - thickness_PML_x

	for ix=1:nx
		# abscissa of current grid point along the damping profile
		xval = δx * real(ix-1)

		#---------- left edge
		if(flags[1]) 

			# define damping profile at the grid points
			abscissa_in_PML = xoriginleft - xval
			if(abscissa_in_PML >= 0.0) 
				abscissa_normalized = abscissa_in_PML / thickness_PML_x
				d_x[ix] = d0_x * abscissa_normalized.^NPOWER
				# this taken from Gedney page 8.2
				k_x[ix] = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized.^NPOWER
				alpha_x[ix] = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.01e0 * ALPHA_MAX_PML
			end

			# define damping profile at half the grid points
			abscissa_in_PML = xoriginleft - (xval + δx/2.e0)
			if(abscissa_in_PML >= 0.0) 
				abscissa_normalized = abscissa_in_PML / thickness_PML_x
				d_x_half[ix] = d0_x * abscissa_normalized.^NPOWER
				# this taken from Gedney page 8.2
				k_x_half[ix] = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized.^NPOWER
				alpha_x_half[ix] = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.01e0 * ALPHA_MAX_PML
			end
		end

		#---------- right edge
		if(flags[2]) 

			# define damping profile at the grid points
			abscissa_in_PML = xval - xoriginright
			if(abscissa_in_PML >= 0.0) 
				abscissa_normalized = abscissa_in_PML / thickness_PML_x
				d_x[ix] = d0_x * abscissa_normalized^NPOWER
				# this taken from Gedney page 8.2
				k_x[ix] = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized^NPOWER
				alpha_x[ix] = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.01e0 * ALPHA_MAX_PML
			end

			# define damping profile at half the grid points
			abscissa_in_PML = xval + δx/2.e0 - xoriginright
			if(abscissa_in_PML >= 0.0) 
				abscissa_normalized = abscissa_in_PML / thickness_PML_x
				d_x_half[ix] = d0_x * abscissa_normalized^NPOWER
				# this taken from Gedney page 8.2
				k_x_half[ix] = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized^NPOWER
				alpha_x_half[ix] = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.01e0 * ALPHA_MAX_PML
			end

		end
#
		# just in case, for -5 at the end
		(alpha_x[ix] < 0.0) ? alpha_x[ix] = 0.0 : nothing
		(alpha_x_half[ix] < 0.0) ? alpha_x_half[ix] = 0.0 : nothing

		# see equation 25 Komatitsch, 2007, get b_x and a_x (need k_x, alpha_x, and d_x before this )
		b_x[ix] = exp(- (d_x[ix] / k_x[ix] + alpha_x[ix]) * δt) 
		b_x_half[ix] = exp(- (d_x_half[ix] / k_x_half[ix] + alpha_x_half[ix]) * δt)

		# this to avoid division by zero outside the PML
		(abs(d_x[ix]) > 1.e-6) ? a_x[ix] = d_x[ix] * (b_x[ix] - 1.e0) / (k_x[ix] * (d_x[ix] + k_x[ix] * alpha_x[ix])) : nothing
		(abs(d_x_half[ix]) > 1.e-6) ? a_x_half[ix] = d_x_half[ix] * (b_x_half[ix] - 1.e0) / (k_x_half[ix] * (d_x_half[ix] + k_x_half[ix] * alpha_x_half[ix])) : nothing

	end
	copyto!(k_xI, inv.(k_x))
	copyto!(k_x_halfI, inv.(k_x_half))



	names=[:a,:b,:kI,:a_half,:b_half,:k_halfI]
	# return [d_x, d_x_half, alpha_x, alpha_x_half]
	return NamedArray([a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI], names)
end



"""
Generate a NamedArray with PML coefficients for all the dimensions that are then stored in the FDTD structs.
"""
function get_pml(mgrid, abs_trbl, args...)
	names=dim_names(length(mgrid))
	return NamedArray([pml_coeff(mgrid[i], [any(abs_trbl .== Symbol(string(dim),"min")), any(abs_trbl .== Symbol(string(dim),"max"))], args...) 
		for (i,dim) in enumerate(names)], names)
end