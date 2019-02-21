mutable struct ParamExptAll
	pa
end


mutable struct PEc
	acqgeom::Acquisition.Geom # mainly positions of receivers are used from here, for the source see P
	datP::Data.TD
	datLP::Data.TD
	tgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}} # tgrid of snapshots
	mgrid::Vector{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}} # model grid of snapshots
end

# type for doing a SE expt
mutable struct ParamExpt
	paTD::Data.P_misfit #  to measure data misfit
	σ::Matrix{Float64} # see Niels et al.
	Q::Matrix{Float64} # Q*k/η (see Niels et al.)
	P::Array{Float64,3} # pressure snapshots (B*P) from seismic modeling
	LP::Array{Float64,3} # ∇⋅([Q*k/η](x,z)∇B*P) 
	ψ::Array{Float64,3} # we are going to record this field
	pad::Param # Param for each snapshot
	pai::Param # Param for each snapshot
end

# allocate 
function ParamExpt(snaps, acqgeom, tgrid, mgrid, Qv, k, η, σ; σobs=nothing, Qobs=nothing)
	length(tgrid)≠size(snaps,3) && error("tgrid and snaps inconsistent")

	P=zeros(size(snaps));
	LP=zeros(size(snaps));
	ψ=zeros(size(snaps));
	for i in eachindex(P)
		P[i]=snaps[i]
	end

	Q= k .* Qv ./ η

	dobs=Data.TD_zeros([:P],tgrid,acqgeom)
	Random.randn!(dobs) # dummy dobs, update later
	paTD=Data.P_misfit(Data.TD_zeros([:P],tgrid,acqgeom),dobs);

	pad=Param(length.(mgrid)...)
	pai=Param(length.(mgrid)...)

	datP=Data.TD_zeros([:P],tgrid,acqgeom)
	datLP=Data.TD_zeros([:P],tgrid,acqgeom)


	pa=ParamExpt(paTD, acqgeom, σ, Q, P, LP, datP, datLP, ψ, pad, pai, tgrid, mgrid)

	# record P and LP
	record!(pa.datP, pa.P, pa.mgrid)
	record!(pa.datLP, pa.LP, pa.mgrid)

	# 
	updateLP!(pa)
	update_observed_data!(pa, σobs,Qobs)
	return pa
end

# given snapshots, record a data obj TD
function record!(data::Data.TD, snaps, mgrid)
	ifield=1;
	iss=1;
	nr=data.acqgeom.nr
	dat=data.d[iss,ifield]
	for ir = 1:nr[iss]
		din=view(dat,:,ir)
		irx=argmin(abs.(mgrid[2].-data.acqgeom.rx[iss][ir]))
		irz=argmin(abs.(mgrid[1].-data.acqgeom.rz[iss][ir]))
		for it in 1:length(data.tgrid)
			din[it]=snaps[irz,irx,it]
		end
	end
end

# generate LP from P, need to optimize this routine, has a loop over snapshots
function updateLP!(pa::ParamExpt)
	updateA!(pa.pad, pa.Q) # put Q in the solver per each snapshot
	@showprogress "updating LP\t" for it in 1:length(pa.tgrid)
		snap_in=view(pa.P,:,:,it)
		snap_out=view(pa.LP,:,:,it)
		solve!(snap_out, snap_in, pa.pad, -1)
	end
	# record LP
	record!(pa.datLP, pa.LP, pa.mgrid)
end

# update psi using LP, need to optimize this routine, has a loop over snapshots
# most expensive routine
function updateψ!(pa::ParamExpt)
	updateA!(pa.pai, pa.σ) # put sigma in the snapshot solver
	qrA=factorize(pa.pai.A)
	@showprogress "updating ψ\t" for it in 1:length(pa.tgrid)
		snap_in=view(pa.LP,:,:,it)
		snap_out=view(pa.ψ,:,:,it)
		solve!(snap_out, snap_in, pa.pai, 1; A=qrA)
	end
end

# do forward modeling to generate observed data in paTD
function update_observed_data!(pa::ParamExpt, σobs, Qobs)
	# save sigma
	σsave=copy(pa.σ)
	Qsave=copy(pa.Q)
	# replace sigma
	copyto!(pa.σ,σobs)
	copyto!(pa.Q,Qobs)
	# do mod
	updateLP!(pa)
	updateψ!(pa)
	# record
	record!(pa.paTD.y, pa.ψ, pa.mgrid)
	# put sigma back
	copyto!(pa.σ, σsave)
	copyto!(pa.Q, Qsave)
end

# do forward modeling, record paTD.x, and measure the misfit in paTD
function funcσ(x, pa::ParamExpt)
	# same as above, but for modelled data obj
	copyto!(pa.σ,x)
	updateψ!(pa)
	record!(pa.paTD.x, pa.ψ, pa.mgrid)
	f = Data.func_grad!(pa.paTD)
	return f
end
# do forward modeling, record paTD.x, and measure the misfit in paTD
function funcQ(x, W, pa::ParamExpt)
	# same as above, but for modelled data obj
	#copyto!(pa.Q,x)
	#for i in eachindex(x)
	#	pa.Q[i]=x[i]
	#end
	#updateLP!(pa)
	#updateψ!(pa)
	#record!(pa.paTD.x, pa.ψ, pa.mgrid)
	#f = Data.func_grad!(pa.paTD)
	y=sum(abs,W*x)
#	f=3.0
	return y
end
