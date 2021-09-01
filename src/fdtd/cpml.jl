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


	d_x=zeros(nx);	d_x_half=zeros(nx)
	k_x=ones(nx);	k_xI=ones(nx)
	k_x_half=ones(nx);	k_x_halfI=ones(nx)
	alpha_x=zeros(nx);	alpha_x_half=zeros(nx)
	a_x=zeros(nx);	a_x_half=zeros(nx)
	b_x=zeros(nx);	b_x_half=zeros(nx)


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



	names=[:a,:b,:kI,:a_half,:b_half,:k_halfI,:d_x,:d_x_half,:alpha_x,:alpha_x_half]
	# return [d_x, d_x_half, alpha_x, alpha_x_half]
	return NamedArray([a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI,d_x, d_x_half, alpha_x,alpha_x_half], names)
end



"""
Generate a NamedArray with PML coefficients for all the dimensions that are then stored in the FDTD structs.
"""
function get_pml(mgrid, abs_trbl, args...)
	names=[:z,:y,:x][1:length(mgrid)]
	pnames=broadcast(x->Symbol("pml",string(x)), names)
	for (i,dim) in enumerate(names) 
		d1=Symbol(string(dim),"min")
		d2=Symbol(string(dim),"max")
		pn=pnames[i]
		@eval $pn=pml_coeff(mgrid[i], [any(abs_trbl .== d1), any(abs_trbl .== d2)], args...)
	end
	# pml_variables
	pml_names=vcat([:a_x,:b_x,:k_xI,:a_x_half,:b_x_half,:k_x_halfI], [:a_z,:b_z,:k_zI,:a_z_half,:b_z_half,:k_z_halfI])
	@eval pml=NamedArray([$(pnames...)], pnames)
	return pml
end