
"""
Return functional and gradient of the LS objective 
* `last_x::Vector{Float64}` : buffer is only updated when x!=last_x, and modified such that last_x=x
"""
function func(x::Vector{Float64}, last_x::Vector{Float64}, pa::Param)
	global to

	if(!isequal(x, last_x))
		copyto!(last_x, x)
		# do forward modelling, apply F x
		@timeit to "F!" F!(pa, x, pa.attrib_mod)
	end

	# compute misfit 
	f = Data.func_grad!(pa.paTD)
	return f
end

function grad!(storage, x::Vector{Float64}, last_x::Vector{Float64}, pa::Param)
	global to

	# (inactive when applied on same model)
	if(!isequal(x, last_x))
		copyto!(last_x, x)
		# do forward modelling, apply F x 
		@timeit to "F!" F!(pa, x, pa.attrib_mod)
	end

	# compute functional and get ∇_d J (adjoint sources)
	f = Data.func_grad!(pa.paTD, :dJx);

	# update adjoint sources after time reversal
	update_adjsrc!(pa.adjsrc, pa.paTD.dJx, pa.adjacqgeom)

	# do adjoint modelling here with adjoint sources Fᵀ F P x
	@timeit to "Fadj!" Fadj!(pa)	

	# adjoint of interpolation
        spray_gradient!(storage,  pa, pa.attrib_mod)

	return storage
end 



function ζfunc(x, last_x, pa::Param, ::LS)
	return func(x, last_x, pa)
end


function ζgrad!(storage, x, last_x, pa::Param, ::LS)
	return grad!(storage, x, last_x, pa)
end


function ζfunc(x, last_x, pa::Param, obj::LS_prior)
	f1=func(x, last_x, pa)

	f2=Misfits.error_squared_euclidean!(nothing, x, pa.mx.prior, pa.mx.w, norm_flag=false)

	return f1*obj.α[1]+f2*obj.α[2]
end

function ζgrad!(storage, x, last_x, pa::Param, obj::LS_prior)
	g1=pa.mx.gm[1]
	grad!(g1, x, last_x, pa)

	g2=pa.mx.gm[2]
	Misfits.error_squared_euclidean!(g2, x, pa.mx.prior, pa.mx.w, norm_flag=false)

	rmul!(g1, obj.α[1])
	rmul!(g2, obj.α[2])
	for i in eachindex(storage)
		@inbounds storage[i]=g1[i]+g2[i]
	end
	return storage
end


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
	copyto!(pa.paTD.x,dcal)
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
	copyto!(pa.paTD.x,dcal)

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


function operator_Born(pa)
	(typeof(pa.attrib_mod) ≠ ModFdtdBorn) && error("need born flag")
	fw=(y,x)->Fborn_map!(y, x, pa)
	bk=(y,x)->Fadj_map!(y, x, pa)

	return LinearMap(fw, bk, 
		  length(pa.paTD.dJx),  # length of output
		  xfwi_ninv(pa), # length of input
		  ismutating=true)
end

function Fborn_map!(y, x, pa)
	F!(pa, x, ModFdtdBorn())
	copyto!(y, pa.paTD.x)
end

function Fadj_map!(y, x, pa)
	copyto!(pa.paTD.dJx, x)

	# adjoint sources
	update_adjsrc!(pa.adjsrc, pa.paTD.dJx, pa.adjacqgeom)

	# adjoint simulation
	Fadj!(pa)

	# adjoint of interpolation
	spray_gradient!(y,  pa, ModFdtdBorn())
end



