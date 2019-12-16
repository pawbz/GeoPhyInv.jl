
"""
Return functional and gradient of the LS objective 
* `last_x::Vector{Float64}` : buffer is only updated when x!=last_x, and mediumified such that last_x=x
"""
function func(x::Vector{Float64}, last_x::Vector{Float64}, pa::PFWI)
	global fwi_to

	if(!isequal(x, last_x))
		copyto!(last_x, x)
		# do forward modelling, apply F x
		@timeit_debug fwi_to "F!" F!(pa, x)
	end

	# compute misfit 
	f = func_grad!(pa.paTD)
	return f
end

function grad!(storage, x::Vector{Float64}, last_x::Vector{Float64}, pa::PFWI)
	global fwi_to
	reset_timer!(fwi_to)

	# (inactive when applied on same model)
	if(!isequal(x, last_x))
		copyto!(last_x, x)
		# do forward modelling, apply F x 
		@timeit_debug fwi_to "F!" F!(pa, x)
	end

	# compute functional and get ∇_d J (adjoint sources)
	f = func_grad!(pa.paTD, :dJx);

	# update adjoint sources after time reversal
	update_adjsrc!(pa.adjsrc, pa.paTD.dJx, pa.adjageom)

	# do adjoint modelling here with adjoint sources Fᵀ F P x
	@timeit_debug fwi_to "Fadj!" Fadj!(pa)	

	# adjoint of interpolation
	spray_gradient!(storage,  pa)

	return storage
end 



function ζfunc(x, last_x, pa::PFWI{T1,T2,T3}) where {T1,T2,T3<:Union{LS,Migr,Migr_FD}}
	return func(x, last_x, pa)
end


function ζgrad!(storage, x, last_x, pa::PFWI{T1,T2,T3})  where {T1,T2,T3<:Union{LS,Migr,Migr_FD}} 
	return grad!(storage, x, last_x, pa)
end


function ζfunc(x, last_x, pa::PFWI{T1,T2,LS_prior})  where {T1,T2} 
	f1=func(x, last_x, pa)

	# calculate the generalized least-squares error
	# note: change the inverse model covariance matrix `pmgls.Q` accordingly
	f2=Misfits.func_grad!(nothing, x, pa.mx.prior, obj.pmgls)

	return f1*obj.pdgls+f2
end

function ζgrad!(storage, x, last_x, pa::PFWI{T1,T2,LS_prior}) where {T1,T2}  
	g1=pa.mx.gm[1]
	grad!(g1, x, last_x, pa)

	g2=pa.mx.gm[2]
	Misfits.func_grad!(g2, x, pa.mx.prior, obj.pmgls)

	rmul!(g1, obj.pdgls)

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
* `pa::PFWI` : parameters that are constant during the inversion 
* if x is absent, using `pa.mediumm` for modeling
"""
function F!(pa::PFWI{Fdtd,T1,T2}, x) where {T1,T2}

	# switch off born scattering
	#pa.paf.c.born_flag=false

	# initialize boundary, as we will record them now
	initialize_boundary!(pa.paf)

	if(!(x===nothing))
		# project x, which lives in mediumi, on to model space (mediumm)
		x_to_mediumm!(pa, x)
	end

	# update model in the forward engine
	update!(pa.paf.c, pa.mediumm)

	pa.paf.c.activepw=[1,]
	pa.paf.c.illum_flag=false
	pa.paf.c.sflags=[2, 0]
	pa.paf.c.rflags=[1, 0] # record only after first scattering
	update_srcwav!(pa.paf,[pa.srcwav,pa.adjsrc])
	pa.paf.c.backprop_flag=1
	pa.paf.c.gmodel_flag=false

	update!(pa.paf);

	# copy data to evaluate misfit
	dcal=pa.paf.c.data[1]
	copyto!(pa.paTD.x,dcal)
end


"""
Born modeling with `mediumm` as the perturbed model and `mediumm0` as the background model.
"""
function F!(pa::PFWI{FdtdBorn,T1,T2}, x) where {T1,T2}

	# update background model in the forward engine 
	update!(pa.paf.c, pa.mediumm0)
	if(!(x===nothing))
		# project x, which lives in mediumi, on to model space (mediumm)
		x_to_mediumm!(pa, x)
	end
	# update perturbed models in the forward engine
	update_δmods!(pa.paf.c, pa.mediumm)

	Fbornmod!(pa::PFWI)
end

function Fbornmod!(pa::PFWI{FdtdBorn,T1,T2}) where {T1,T2} 

	# switch on born scattering
	#pa.paf.c.born_flag=true

	pa.paf.c.activepw=[1,2] # two wavefields are active
	pa.paf.c.illum_flag=false 
	pa.paf.c.sflags=[2, 0] # no sources on second wavefield
	pa.paf.c.rflags=[0, 1] # record only after first scattering

	# source wavelets (for second wavefield, they are dummy)
	update_srcwav!(pa.paf,[pa.srcwav,pa.adjsrc])

	# actually, should record only when background field is changed
	pa.paf.c.backprop_flag=1 # store boundary values for gradient later

	pa.paf.c.gmodel_flag=false # no gradient

	update!(pa.paf);
	dcal=pa.paf.c.data[2]
	copyto!(pa.paTD.x,dcal)

	# switch off born scattering once done
	#pa.paf.c.born_flag=false
end


"""
Perform adjoint modelling in `paf` using adjoint sources `adjsrc`.
"""
function Fadj!(pa::PFWI)

	# need to explicitly turn off the born flag for adjoint modelling
	#pa.paf.c.born_flag=false

	# require gradient, switch on the flag
	pa.paf.c.gmodel_flag=true

	# both wavefields are active
	pa.paf.c.activepw=[1,2]

	# no need of illum during adjoint modeling
	pa.paf.c.illum_flag=false

	# force boundaries in first pw and back propagation for second pw
	pa.paf.c.sflags=[-2,2] 
	pa.paf.c.backprop_flag=-1

	# update source wavelets in paf using adjoint sources
	update_srcwav!(pa.paf,[pa.srcwav,pa.adjsrc])

	# no need to record data during adjoint propagation
	pa.paf.c.rflags=[0,0]

	# adjoint modelling
	update!(pa.paf);

	# put rflags back
	pa.paf.c.rflags=[1,1]

	return pa.paf.c.gradient
end


"""
```
F=LinearMap(pa)
```
If `pa` is an instance of `SeisInvExpt`, then 
return the linearized forward modeling operator `F`, such that
`F*x` can be computed without explicitly storing the operator matrix (see `LinearMaps.jl`).
The imaging/migration operator is given by `transpose(F)`. 
These operators are the building blocks of iterative optimization schemes.
"""
function LinearMaps.LinearMap(pa::PFWI{FdtdBorn,T2,T3}) where {T2,T3}
	fw=(y,x)->Fborn_map!(y, x, pa)
	bk=(y,x)->Fadj_map!(y, x, pa)

	return LinearMap(fw, bk, 
		sum(length.(pa.paTD.dJx)),  # length of output
		  xfwi_ninv(pa), # length of input
		  ismutating=true)
end

function Fborn_map!(δy, δx, pa)
	δx_to_δmods!(pa, δx)
	Fbornmod!(pa)
	copyto!(δy, pa.paTD.x)
end

function Fadj_map!(δy, δx, pa)
	copyto!(pa.paTD.dJx, δx)

	# adjoint sources
	update_adjsrc!(pa.adjsrc, pa.paTD.dJx, pa.adjageom)

	# adjoint simulation
	Fadj!(pa)

	# chain rule corresponding to reparameterization
	pert_chainrule!(pa.mxm.gx, pa.paf.c.gradient, pa.mediumm0, pa.parameterization)

	# finally, adjoint of interpolation
	Interpolation.interp_spray!(δy, 
			     pa.mxm.gx, pa.paminterp, :spray, 
			     count(pa.parameterization.≠:null))
end



