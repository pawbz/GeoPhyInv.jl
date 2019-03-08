struct Fσ end
struct FGσ end

# type for doing a SE expt
mutable struct ParamExpt{T}
#	paTD::Data.P_misfit #  to measure data misfit
#	acqgeom::Acquisition.Geom # mainly positions of receivers are used from here, for the source see P
	σ::Matrix{T} # see Niels et al.
	Q::Matrix{T} # Q*k/η (see Niels et al.)
	P::Array{T,3} # pressure snapshots (B*P) from seismic modeling
	LP::Array{T,3} # ∇⋅([Q*k/η](x,z)∇B*P) 
	data::Array{T,2}
	data_obs::Array{T,2}
	data_misfit::Array{T,1}
#	datP::Data.TD
#	datLP::Data.TD
	ψ::Vector{T} # state variable, we are going to record this state data=ACQ*ψ  
	adjsrc::Vector{T}
	adjψ::Vector{T} # adjoint state variable
	paQ::Param{T} # Param for each snapshot
	paσ::Param{T} # Param for each snapshot
	tgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}} # tgrid of snapshots
	mgrid::Vector{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}} # model grid of snapshots
	ACQ::SparseMatrixCSC{T,Int64}
	g::Matrix{T}
	gtemp::Matrix{T}
end

# allocate 
function ParamExpt(snaps, tgrid, mgrid, Qv, k, η, σ; σobs=nothing, Qobs=nothing)
	nz=length(mgrid[1])
	nx=length(mgrid[2])
	nt=length(tgrid)
	(size(snaps)≠(nz,nx,nt)) && error("tgrid, mgrid and snaps inconsistent")
	(size(Qv)≠(nz,nx)) && error("tgrid, mgrid and Qv inconsistent")
	(size(k)≠(nz,nx)) && error("tgrid, mgrid and k inconsistent")
	(size(η)≠(nz,nx)) && error("tgrid, mgrid and η inconsistent")
	(size(σ)≠(nz,nx)) && error("tgrid, mgrid and σ inconsistent")

	LP=zero(snaps);
	ψ=zeros(nz*nx);	
	adjsrc=zeros(nz*nx);	
	adjψ=zeros(nz*nx+1);	

	Q= k .* Qv ./ η # combine all these

#	dobs=Data.TD_zeros([:P],tgrid,acqgeom)
	#Random.randn!(dobs) # dummy dobs, update later
#	paTD=Data.P_misfit(Data.TD_zeros([:P],tgrid,acqgeom),dobs);

	paQ=Param(mgrid, Q)
	paσ=Param(mgrid, σ)

#	datP=Data.TD_zeros([:P],tgrid,acqgeom)
#	datLP=Data.TD_zeros([:P],tgrid,acqgeom)


# 5 receivers
nr=nz*nx
	ACQ=spzeros(nr,nz*nx)
	for i in 1:2:nr*nz*nx
		ACQ[i]=randn()
	end
	data=zeros(nt,size(ACQ,1))
	data_obs=randn(nt,size(ACQ,1))
	data_misfit=randn(size(ACQ,1))

	pa=ParamExpt(#paTD, 
	      #acqgeom, 
	      σ, Q, snaps,  LP, #
	      data,
	      data_obs,
	      data_misfit,
	      ψ,
	      adjsrc,
	      adjψ,
	     # ssnaps,
	     # data,
	     #datP, datLP, 
	     #ψ, 
	     paQ, paσ, 
	     tgrid, mgrid,
	     ACQ,
	     zero(σ),zero(σ))

	# record P and LP
#	record!(pa.datP, pa.P, pa.mgrid)
#	record!(pa.datLP, pa.LP, pa.mgrid)

	# 
	if(!(σobs===nothing) && !(Qobs===nothing))
		updateLP!(pa, Qobs)
		update_observed_data!(pa, σobs)
	end
	return pa
end

#=
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
=#



# generate LP from P, need to optimize this routine, has a loop over snapshots
function updateLP!(pa::ParamExpt, Q)
	updateA!(pa.paQ, Q)
	@showprogress "updating LP\t" for it in 1:length(pa.tgrid)
		snap_in=view(pa.P,:,:,it)
		snap_out=view(pa.LP,:,:,it)
		applyA!(snap_out, snap_in, pa.paQ; A=pa.paQ.A)
	end
	return pa
end

# forward modeling to generate ψ
# update psi using LP, has a loop over snapshots
# most expensive routine
function mod!(pa::ParamExpt, σ, mode)
	updateA!(pa.paσ, σ)
	fill!(pa.g, 0.0)
	@showprogress "updating ψ\t" for it in 1:length(pa.tgrid)
		snap_in=view(pa.LP,:,:,it)
		forwψ!(pa.ψ,snap_in,pa.paσ)

		# record
		dat_slice=view(pa.data,it,:)
		mul!(dat_slice,pa.ACQ,pa.ψ) # record

		# wanna do adjoint modeling? depends on mode
		dat_slice_obs=view(pa.data_obs,it,:)
		adjψ!(pa.adjψ, pa.ψ, pa.adjsrc, pa.g, pa.ACQ, pa.data_misfit, dat_slice_obs, pa.paσ, mode)
	end
	return pa
end

function forwψ!(ψ, src, pa::Param)
	applyinvA!(ψ, src, pa; A=pa.A)
end
function adjψ!(adjψ, ψ, adjsrc, g, ACQ, data_misfit, data_obs, pa::Param, ::Fσ)
end
function adjψ!(adjψ, ψ, adjsrc, g, ACQ, data_misfit, data_obs, pa::Param, ::FGσ)

	mul!(data_misfit,ACQ,ψ)
	for i in eachindex(data_misfit)
		data_misfit[i] -= data_obs[i]
	end
	mul!(adjsrc, transpose(ACQ), data_misfit)
	rmul!(adjsrc,-2.0)

	applyinvAt!(adjψ, adjsrc, pa)

	gtemp=vec((adjψ*ψ'))'*(pa.dAdx)

	for i in eachindex(gtemp)
                       g[i]+=gtemp[i]
                 end
end

function func(x,pa::ParamExpt)
	mod!(pa, x, Fσ())
	f=sum(abs2,pa.data.-pa.data_obs)
	return f
end

function grad(x,pa::ParamExpt)
	mod!(pa, x, FGσ())
	return paE.g
end

# do forward modeling to generate observed data in paTD
function update_observed_data!(pa::ParamExpt, σobs)
	# save sigma
	σsave=copy(pa.σ)
	mod!(pa, σobs, Fσ())
	copyto!(pa.data_obs, pa.data)
	# put sigma back
	copyto!(pa.σ, σsave)
	# dummy modeling, just to remove traces of σobs
	mod!(pa, σobs, Fσ())
end

