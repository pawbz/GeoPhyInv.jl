

@inbounds @fastmath function boundary_force_snap_p!(issp::Int64,pac::Paramc,pass::Vector{Paramss},pap::Paramp)
	p=pap.p
	boundary=pass[issp].boundary[5]
	ps=view(p[1],:,:,1)
	bs=view(boundary,:,:,1)
	copy!(ps,bs)
end
@inbounds @fastmath function boundary_force_snap_vxvz!(issp::Int64,pac::Paramc,pass::Vector{Paramss},pap::Paramp)
	# initial conditions from boundary for first propagating field only
	p=pap.p
	pw=p[1]
	boundary=pass[issp].boundary[5]
	ps=view(pw,:,:,2)
	bs=view(boundary,:,:,2)
	copy!(ps,bs)
	ps=view(pw,:,:,3)
	bs=view(boundary,:,:,3)
	copy!(ps,bs)
end
@fastmath @inbounds function boundary_force!(it::Int64,issp::Int64,pac::Paramc,pass::Vector{Paramss},pap::Paramp)
	boundary=pass[issp].boundary
	p=pap.p
	ibx0=pac.ibx0; ibz0=pac.ibz0; ibx1=pac.ibx1; ibz1=pac.ibz1
	boundaryf_l!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	boundaryf_r!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	boundaryf_t!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	boundaryf_b!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
end
@fastmath @inbounds function boundaryf_l!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	pw=p[1]
	for ix=1:3
		@simd for iz=1:ibz1-ibz0+1
			pw[ibz0+iz-1,ibx0+ix-1,1] = boundary[4][iz,ix,it]
		end
	end
end
@fastmath @inbounds function boundaryf_r!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	pw=p[1]
	for ix=1:3
		@simd for iz=1:ibz1-ibz0+1
			pw[ibz0+iz-1,ibx1-ix+1,1] = boundary[2][iz,ix,it]
		end
	end
end
@fastmath @inbounds function boundaryf_t!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	pw=p[1]
	@simd for ix=1:ibx1-ibx0+1
		for iz=1:3
			pw[ibz0+iz-1,ibx0+ix-1,1] = boundary[1][iz,ix,it]
		end
	end
end
@fastmath @inbounds function boundaryf_b!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	pw=p[1]
	@simd for ix=1:ibx1-ibx0+1
		for iz=1:3
			pw[ibz1-iz+1,ibx0+ix-1,1] = boundary[3][iz,ix,it]
		end
	end
end


@inbounds @fastmath function boundary_save_snap_p!(issp::Int64,pac::Paramc,pass::Vector{Paramss},pap::Paramp)
	p=pap.p
	pw=p[1]
	boundary=pass[issp].boundary[5]
	ps=view(pw,:,:,1)
	bs=view(boundary,:,:,1)
	copy!(bs, ps)
end
@inbounds @fastmath function boundary_save_snap_vxvz!(issp::Int64,pac::Paramc,pass::Vector{Paramss},pap::Paramp)
	p=pap.p
	pw=p[1]
	boundary=pass[issp].boundary[5]
	#vx
	bs=view(boundary,:,:,2)
	ps=view(pw,:,:,2)
	copy!(bs, ps)
	scale!(bs,-1.)
	# vz
	bs=view(boundary,:,:,3)
	ps=view(pw,:,:,3)
	copy!(bs, ps)
	scale!(bs,-1.)
end
@fastmath @inbounds function boundary_save!(it::Int64,issp::Int64,pac::Paramc,pass::Vector{Paramss},pap::Paramp)
	boundary=pass[issp].boundary
	p=pap.p
	ibx0=pac.ibx0; ibz0=pac.ibz0; ibx1=pac.ibx1; ibz1=pac.ibz1
	boundarys_l!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	boundarys_r!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	boundarys_t!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	boundarys_b!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
end
@fastmath @inbounds function boundarys_l!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	pw=p[1]
	for ix=1:3
		@simd for iz=1:ibz1-ibz0+1
			boundary[4][iz,ix,it] = pw[ibz0+iz-1,ibx0+ix-1,1] 
		end
	end
end
@fastmath @inbounds function boundarys_r!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	pw=p[1]
	for ix=1:3
		@simd for iz=1:ibz1-ibz0+1
			boundary[2][iz,ix,it] = pw[ibz0+iz-1,ibx1-ix+1,1]
		end
	end
end
@fastmath @inbounds function boundarys_t!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	pw=p[1]
	@simd for ix=1:ibx1-ibx0+1
		for iz=1:3
			boundary[1][iz,ix,it] = pw[ibz0+iz-1,ibx0+ix-1,1]
		end
	end
end
@fastmath @inbounds function boundarys_b!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	pw=p[1]
	@simd for ix=1:ibx1-ibx0+1
		for iz=1:3
			boundary[3][iz,ix,it] = pw[ibz1-iz+1,ibx0+ix-1,1]
		end
	end
end


