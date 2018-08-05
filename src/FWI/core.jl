

"""
Perform a forward simulation.
Update `pa.paTD.x`. 
This simulation is common for both functional and gradient calculation.
During the computation of the gradient, we need an adjoint simulation.
Update the buffer, which consists of the modelled data
and boundary values for adjoint calculation.

# Arguments

* `x::Vector{Float64}` : inversion variable
* `pa::Param` : parameters that are constant during the inversion 
* if x is absent, using `pa.modm` for modeling
"""
function F!(pa::Param, x, ::ModFdtd)

	# initialize boundary, as we will record them now
	Fdtd.initialize_boundary!(pa.paf)

	if(!(x===nothing))
		# project x, which lives in modi, on to model space (modm)
		x_to_modm!(pa, x)
	end

	# update model in the forward engine
	Fdtd.update_model!(pa.paf.c, pa.modm)

	pa.paf.c.activepw=[1,]
	pa.paf.c.illum_flag=false
	pa.paf.c.sflags=[2, 0]
	Fdtd.update_acqsrc!(pa.paf,[pa.acqsrc,pa.adjsrc])
	pa.paf.c.backprop_flag=1
	pa.paf.c.gmodel_flag=false

	Fdtd.mod!(pa.paf);

	# copy data to evaluate misfit
	dcal=pa.paf.c.data[1]
	copy!(pa.paTD.x,dcal)
end


"""
Born modeling with `modm` as the perturbed model and `modm0` as the background model.
"""
function F!(pa::Param, x, ::ModFdtdBorn)
	# switch on born scattering
	pa.paf.c.born_flag=true

	if(!(x===nothing))
		# project x, which lives in modi, on to model space (modm)
		x_to_modm!(pa, x)
	end

	# update background and perturbed models in the forward engine
	Fdtd.update_model!(pa.paf.c, pa.modm0, pa.modm)

	pa.paf.c.activepw=[1,2] # two wavefields are active
	pa.paf.c.illum_flag=false 
	pa.paf.c.sflags=[2, 0] # no sources on second wavefield
	pa.paf.c.rflags=[0, 1] # record only after first scattering

	# source wavelets (for second wavefield, they are dummy)
	Fdtd.update_acqsrc!(pa.paf,[pa.acqsrc,pa.adjsrc])

	pa.paf.c.backprop_flag=1 # store boundary values for gradient later
	pa.paf.c.gmodel_flag=false # no gradient
	

	Fdtd.mod!(pa.paf);
	dcal=pa.paf.c.data[2]
	copy!(pa.paTD.x,dcal)

	pa.paf.c.born_flag=false
end


"""
Perform adjoint modelling in `paf` using adjoint sources `adjsrc`.
"""
function Fadj!(pa::Param)

	# require gradient
	pa.paf.c.gmodel_flag=true

	# both wavefields are active
	pa.paf.c.activepw=[1,2]

	# no need of illum during adjoint modeling
	pa.paf.c.illum_flag=false

	# need to explicitly turn off the born flag for adjoint modelling
	pa.paf.c.born_flag=false

	# force boundaries in first pw and back propagation for second pw
	pa.paf.c.sflags=[3,2] 
	pa.paf.c.backprop_flag=-1

	# update source wavelets in paf using adjoint sources
	Fdtd.update_acqsrc!(pa.paf,[pa.acqsrc,pa.adjsrc])


	# no need to record data
	pa.paf.c.rflags=[0,0]

	# adjoint modelling
	Fdtd.mod!(pa.paf);

	# put rflags back
	pa.paf.c.rflags=[1,1]

	return pa.paf.c.gradient
end

"""
xx = Fadj * Fborn * x
"""
function Fadj_Fborn_x!(xx, x, pa)

	# put  x into pa
	x_to_modm!(pa, x)

	# generate Born Data
	(pa.paf.c.born_flag==false) && error("need born flag")
	Fborn!(pa)

	#f = Data.func_grad!(pa.paTD, :dJx)

	# adjoint sources
	update_adjsrc!(pa.adjsrc, pa.paTD.x, pa.adjacqgeom)

	# adjoint simulation
	Fadj!(pa)

	# adjoint of interpolation
	Seismic_gx!(pa.gmodm, pa.modm, pa.gmodi, pa.modi, xx, pa, 1) 

	return pa
end


function Fadj_Fborn(pa)
	(pa.paf.c.born_flag==false) && error("need born flag")
	fw=(y,x)->F_kernel!(y, x, pa)
	bk=(y,x)->Fadj_kernel!(y, x, pa)

	return LinearMap(fw, bk, 
		  length(pa.paTD.dJx),  # length of output
		  xfwi_ninv(pa), # length of input
		  ismutating=true)
end

function F_kernel!(y, x, pa)
	F!(pa, x, pa.attrib_mod)
	copy!(y, pa.paTD.x)
end

function Fadj_kernel!(y, x, pa)
	copy!(pa.paTD.dJx, x)

	# adjoint sources
	update_adjsrc!(pa.adjsrc, pa.paTD.dJx, pa.adjacqgeom)

	# adjoint simulation
	Fadj!(pa)

	# adjoint of interpolation
        gmodm_to_gx!(y, pa.gmodm, pa, pa.attrib_mod)
end



