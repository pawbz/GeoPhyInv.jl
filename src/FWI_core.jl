

"""
Perform a forward simulation.
This simulation is common for both functional and gradient calculation.
During the computation of the gradient, we need an adjoint simulation.
Update the buffer, which consists of the modelled data
and boundary values for adjoint calculation.

# Arguments

* `x::Vector{Float64}` : inversion variable
* `last_x::Vector{Float64}` : buffer is only updated when x!=last_x, and modified such that last_x=x
* `pa::Param` : parameters that are constant during the inversion 
* `modm::Models.Seismic` : 
"""
function F!(pa::Param, x, last_x=[0.0])
	if(x!=last_x)
		copy!(last_x, x)

		if(!(x===nothing))
			# project x, which lives in modi, on to model space (modm)
			Seismic_x!(pa.modm, pa.modi, x, pa, -1)		

			# update model in the forward engine
			Fdtd.update_model!(pa.paf.c, pa.modm)
		end

		pa.paf.c.activepw=[1,]
		pa.paf.c.illum_flag=false
		pa.paf.c.sflags=[2, 0]
		Fdtd.update_acqsrc!(pa.paf,[pa.acqsrc,pa.adjsrc])
		pa.paf.c.backprop_flag=1
		pa.paf.c.gmodel_flag=false

		Fdtd.mod!(pa.paf);
		dcal=pa.paf.c.data[1]
		copy!(pa.paTD.x,dcal)
	end
end

"""
Born modeling with `modm` as the perturbed model and `modm0` as the background model.
"""
function Fborn!(pa::Param)
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
	
	# switch on born scattering
	pa.paf.c.born_flag=true

	Fdtd.mod!(pa.paf);
	dcal=pa.paf.c.data[2]
	copy!(pa.paTD.x,dcal)

end


"""
Perform adjoint modelling in `paf` using adjoint sources `adjsrc`.
"""
function Fadj!(pa::Param)

	# both wavefields are active
	pa.paf.c.activepw=[1,2]

	# no need of illum during adjoint modeling
	pa.paf.c.illum_flag=false

	# force boundaries in first pw and back propagation for second pw
	pa.paf.c.sflags=[3,2] 
	pa.paf.c.backprop_flag=-1

	# update source wavelets in paf using adjoint sources
	Fdtd.update_acqsrc!(pa.paf,[pa.acqsrc,pa.adjsrc])

	# require gradient
	pa.paf.c.gmodel_flag=true

	# adjoint modelling
	Fdtd.mod!(pa.paf);

	# copy gradient
	copy!(pa.gmodm,pa.paf.c.gmodel)
end

"""
xx = Fadj * Fborn * x
"""
function Fadj_Fborn_x!(xx, x, pa)

	# put  x into pa
	Seismic_x!(pa.modm, pa.modi, x, pa, -1)		

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


