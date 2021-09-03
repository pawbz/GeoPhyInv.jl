
function get_boundary_indices(mgrid,::FdtdOld)
	# where to store the boundary values (careful, born scaterrers cannot be inside these!?)
	ibx0=npml-2; ibx1=length(mgrid[2])+npml+3
	ibz0=npml-2; ibz1=length(mgrid[1])+npml+3
#	ibx0=npml+1; ibx1=length(mgrid[2])+npml
#	ibz0=npml+1; ibz1=length(mgrid[1])+npml
#	println("**** Boundary Storage Changed **** ")

	# for snaps
	isx0, isz0=npml, npml

	return NamedArray([ibx0,ibx1,ibz0,ibz1,isx0,isz0],([:bx0,:bx1,:bz0,:bz1,:sx0,:sz0],))
	
end

function get_boundary_indices(mgrid,::FdtdElastic)
	return NamedArray([1],[:a]) # some dummy, update later	
end



@inbounds @fastmath function boundary_force_snap_p!(issp::Int64,pac,pap)
	ps=pap[1].w1[:t][:p]
	boundary=pap[1].ss[issp].boundary[5]
	bs=view(boundary,:,:,1)
	copyto!(ps,bs)
end
@inbounds @fastmath function boundary_force_snap_vxvz!(issp::Int64,pac,pap)
	# initial conditions from boundary for first propagating field only
	ps=pap[1].w1[:t][:vx]
	boundary=pap[1].ss[issp].boundary[5]
	bs=view(boundary,:,:,2)
	copyto!(ps,bs)
	ps=pap[1].w1[:t][:vz]
	bs=view(boundary,:,:,3)
	copyto!(ps,bs)
end
@fastmath @inbounds function boundary_force!(it::Int64,issp::Int64,pac,pap)
	boundary=pap[1].ss[issp].boundary
	p=pap[1].w1[:t]
	ibx0=pac.bindices[:bx0]; ibz0=pac.bindices[:bz0]; ibx1=pac.bindices[:bx1]; ibz1=pac.bindices[:bz1]
	boundaryf_l!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	boundaryf_r!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	boundaryf_t!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	boundaryf_b!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
end
@fastmath @inbounds function boundaryf_l!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	pw=p[:p]
	for ix=1:3
		@simd for iz=1:ibz1-ibz0+1
			pw[ibz0+iz-1,ibx0+ix-1] = boundary[4][iz,ix,it]
		end
	end
end
@fastmath @inbounds function boundaryf_r!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	pw=p[:p]
	for ix=1:3
		@simd for iz=1:ibz1-ibz0+1
			pw[ibz0+iz-1,ibx1-ix+1] = boundary[2][iz,ix,it]
		end
	end
end
@fastmath @inbounds function boundaryf_t!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	pw=p[:p]
	@simd for ix=1:ibx1-ibx0+1
		for iz=1:3
			pw[ibz0+iz-1,ibx0+ix-1] = boundary[1][iz,ix,it]
		end
	end
end
@fastmath @inbounds function boundaryf_b!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	pw=p[:p]
	@simd for ix=1:ibx1-ibx0+1
		for iz=1:3
			pw[ibz1-iz+1,ibx0+ix-1] = boundary[3][iz,ix,it]
		end
	end
end


@inbounds @fastmath function boundary_save_snap_p!(issp::Int64,pac,pap)
	ps=pap[1].w1[:t][:p]
	boundary=pap[1].ss[issp].boundary[5]
	bs=view(boundary,:,:,1)
	copyto!(bs, ps)
end
@inbounds @fastmath function boundary_save_snap_vxvz!(issp::Int64,pac,pap)
	boundary=pap[1].ss[issp].boundary[5]
	#vx
	bs=view(boundary,:,:,2)
	ps=pap[1].w1[:t][:vx]
	copyto!(bs, ps)
	rmul!(bs,-1.)
	# vz
	bs=view(boundary,:,:,3)
	ps=pap[1].w1[:t][:vz]
	copyto!(bs, ps)
	rmul!(bs,-1.)
end
@fastmath @inbounds function boundary_save!(it::Int64,issp::Int64,pac,pap)
	boundary=pap[1].ss[issp].boundary
	p=pap[1].w1[:t]
	ibx0=pac.bindices[:bx0]; ibz0=pac.bindices[:bz0]; ibx1=pac.bindices[:bx1]; ibz1=pac.bindices[:bz1]
	boundarys_l!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	boundarys_r!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	boundarys_t!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	boundarys_b!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
end
@fastmath @inbounds function boundarys_l!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	pw=p[:p]
	for ix=1:3
		@simd for iz=1:ibz1-ibz0+1
			boundary[4][iz,ix,it] = pw[ibz0+iz-1,ibx0+ix-1] 
		end
	end
end
@fastmath @inbounds function boundarys_r!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	pw=p[:p]
	for ix=1:3
		@simd for iz=1:ibz1-ibz0+1
			boundary[2][iz,ix,it] = pw[ibz0+iz-1,ibx1-ix+1]
		end
	end
end
@fastmath @inbounds function boundarys_t!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	pw=p[:p]
	@simd for ix=1:ibx1-ibx0+1
		for iz=1:3
			boundary[1][iz,ix,it] = pw[ibz0+iz-1,ibx0+ix-1]
		end
	end
end
@fastmath @inbounds function boundarys_b!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	pw=p[:p]
	@simd for ix=1:ibx1-ibx0+1
		for iz=1:3
			boundary[3][iz,ix,it] = pw[ibz1-iz+1,ibx0+ix-1]
		end
	end
end


