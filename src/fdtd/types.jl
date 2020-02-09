
# this is extensively used to stack arrays
# define a specific namedarray
NamedStack{T}=NamedArray{T,1,Array{T,1},Tuple{OrderedCollections.OrderedDict{Symbol,Int64}}}




"""
Modelling parameters per every supersource for each worker
"""
mutable struct P_x_worker_x_pw_x_ss
	iss::Int64
	wavelets::Vector{NamedStack{Vector{Float64}}} # [ .. for it in 1:nt]
	ssprayw::Matrix{Float64}
	records::NamedStack{Matrix{Float64}}
	rinterpolatew::Matrix{Float64}
	sindices::NamedStack{Vector{Int64}} # contains [:x1,:x2,:z1,:z2]
	rindices::NamedStack{Vector{Int64}}
    boundary::Vector{Array{Float64,3}}
	snaps::Array{Float64,3}
	illum::Matrix{Float64}
	grad_mod::NamedStack{Matrix{Float64}} # w.r.t different coeffs
	#=
	grad_modtt::Matrix{Float64} 
	grad_modrrvx::Matrix{Float64}
	grad_modrrvz::Matrix{Float64}
	=#
end


"""
PFdtdeters per every worker, not necessarily for every supersource.
Note that a single worker can take care of multiple supersources.
"""
mutable struct P_x_worker_x_pw
	ss::Vector{P_x_worker_x_pw_x_ss}
	w2::NamedStack{NamedStack{Matrix{Float64}}} # p, vx, vz
	#=
	pp::NamedStack{Matrix{Float64}} # same as above, at previous time step
	ppp::NamedStack{Matrix{Float64}} # ../../p, vx, vz
	dpdx::NamedStack{Matrix{Float64}} # x derivatives of p, vx, vz
	dpdz::NamedStack{Matrix{Float64}} # z derivatives of p, vx, vz
	=#
	memory_pml::NamedStack{Matrix{Float64}} # PML related (stored dpdx, dpdz, dvxdx, dvzdz)
	born_svalue_stack::Matrix{Float64} # used for born modeling 
end


"""
Modelling parameters common for all supersources
# Keyword Arguments that are modified by the method (some of them are returned as well)

* `gradient::Vector{Float64}=Medium(model.mgrid)` : gradient model modified only if `gmodel_flag`
* `TDout::Vector{Data.TD}=[Data.TD_zeros(rfields,tgridmod,ageom[ip]) for ip in 1:length(findn(rflags))]`
* `illum::Array{Float64,2}=zeros(length(model.mgrid[1]), length(model.mgrid[2]))` : source energy if `illum_flag`
* `boundary::Array{Array{Float64,4},1}` : stored boundary values for first propagating wavefield 
* `snaps::Array{Float64,4}=zeros(length(model.mgrid[1]),length(model.mgrid[2]),length(tsnaps),ageom[1].nss)` :snapshots saved at `tsnaps`

# Return (in order)

* modelled data for each propagating wavefield as `Vector{TD}`
* stored boundary values of the first propagating wavefield as `Array{Array{Float64,4},1}` (use for backpropagation)
* final conditions of the first propagating wavefield as `Array{Float64,4}` (use for back propagation)
* gradient model as `Seismic`
* stored snaps shots at tsnaps as Array{Float64,4} 
"""
mutable struct P_common{T}
	jobname::Symbol
	attrib_mod::T
	activepw::Vector{Int64}
	exmodel::Medium
	model::Medium
	ageom::Vector{AGeom}
	srcwav::Vector{SrcWav}
	abs_trbl::Vector{Symbol}
	sfields::Vector{Vector{Symbol}}
	isfields::Vector{Vector{Int64}}
	sflags::Vector{Int64} 
	rfields::Vector{Symbol}
	rflags::Vector{Int64}
	fc::NamedStack{Float64}
	ic::NamedStack{Int64}
	pml::NamedStack{Vector{Float64}}
	#=
	a_x::Vector{Float64}
	b_x::Vector{Float64}
	k_xI::Vector{Float64}
	a_x_half::Vector{Float64}
	b_x_half::Vector{Float64}
	k_x_halfI::Vector{Float64} 
	a_z::Vector{Float64}
	b_z::Vector{Float64}
	k_zI::Vector{Float64}
	a_z_half::Vector{Float64}
	b_z_half::Vector{Float64}
	k_z_halfI::Vector{Float64}
	=#
	mod::NamedStack{Matrix{Float64}}
	mod3d::NamedStack{Array{Float64,3}}
	#=
	modtt::Matrix{Float64}
	modttI::Matrix{Float64} # just storing inv(modtt) for speed
	modrr::Matrix{Float64}
	modrrvx::Matrix{Float64}
	modrrvz::Matrix{Float64}
	=#
	δmod::NamedStack{Matrix{Float64}}
	δmodall::Vector{Float64}
	#=
	δmodtt::Matrix{Float64}
	δmodrr::Matrix{Float64} 
	δmodrrvx::Matrix{Float64}
	δmodrrvz::Matrix{Float64}
	δmod::Vector{Float64} # perturbation vector (KI, rhoI)
	=#
	gradient::Vector{Float64}  # output gradient vector w.r.t (KI, rhoI)
	grad_mod::NamedStack{SharedArrays.SharedArray{Float64,2}}
	#=
	grad_modtt_stack::SharedArrays.SharedArray{Float64,2} # contains gmodtt
	grad_modrrvx_stack::SharedArrays.SharedArray{Float64,2}
	grad_modrrvz_stack::SharedArrays.SharedArray{Float64,2}
	grad_modrr_stack::Array{Float64,2}
	=#
	illum_flag::Bool
	illum_stack::SharedArrays.SharedArray{Float64,2}
	backprop_flag::Int64
	snaps_flag::Bool
	itsnaps::Vector{Int64}
	gmodel_flag::Bool
	bindices::NamedStack{Int64}
	#=
	ibx0::Int64
	ibz0::Int64
	ibx1::Int64
	ibz1::Int64
	isx0::Int64
	isz0::Int64
	=#
	datamat::SharedArrays.SharedArray{Float64,3}
	data::Vector{Data}
	# attenuation related parameters
	verbose::Bool
end 

