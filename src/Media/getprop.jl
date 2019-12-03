"""
Get other dependent model parameters of a seismic model
that are not present in `Seismic`.

* `:rhoI` : inverse of density
* `:Zp` : P-wave impedance
"""
function Base.getindex(mod::Medium, s::Symbol)
	nznx=length(mod.m[:vp])
	vp=mod.m[:vp]; rho=mod.m[:rho]
	(s ∉ names(mod.m)[1]) && (mod.m=combine(mod.m,NamedArray([zero(vp)],[[s],])))
	x=mod.m[s]
	if(s == :χrhoI)
		@inbounds for i in 1:nznx; x[i]=χ(inv(rho[i]),mod.ref[:rhoI],1); end
	elseif(s == :χrho)
		@inbounds for i in 1:nznx; x[i]=χ(rho[i],mod.ref[:rho],1); end
	elseif(s == :χvp)
		@inbounds for i in 1:nznx; x[i]=χ(vp[i],mod.ref[:vp],1); end
	elseif(s == :χK)
		@inbounds for i in 1:nznx; x[i]=χ(vp[i]*vp[i]*rho[i],mod.ref[:K],1); end
	elseif(s == :χmu)
		@inbounds for i in 1:nznx; x[i]=χ(mod.m[:vs][i]*mod.m[:vs][i]*rho[i],
			mod.ref[:mu],1); end
	elseif(s == :χKI)
		@inbounds for i in 1:nznx; x[i]=χ(inv(vp[i]*vp[i]*rho[i]),mod.ref[:KI],1); end
	elseif(s == :KI)
		@inbounds for i in 1:nznx; x[i]=inv(vp[i]*vp[i]*rho[i]); end
	elseif(s == :K)
		@inbounds for i in 1:nznx; x[i]=vp[i]*vp[i]*rho[i]; end
	elseif(s == :Zp)
		@inbounds for i in 1:nznx; x[i]=vp[i]*rho[i]; end
	end
	return mod.m[s]
end



#=



function Seismic_get(mod::Seismic, attrib::Symbol)
	# allocate
	rout=zeros(length(mod.mgrid[1]), length(mod.mgrid[2]))
	Seismic_get!(rout, mod, [attrib])
	return rout
end
function Seismic_get!(x, mod::Seismic, attribvec::Vector{Symbol})
	rho=mod.χrho; χ!(rho, mod.ref.rho,-1) # undo it later
	vp=mod.χvp; χ!(vp, mod.ref.vp,-1) # undo it later
	vs=mod.χvs; χ!(vs, mod.ref.vs,-1) # undo it later
	nznx=length(rho)
	(length(x)≠(count(attribvec.≠ :null)*nznx)) &&  error("size x")
	i0=0; 
	for attrib in attribvec
		if(attrib ≠:null)
		if(attrib == :rhoI)
			@inbounds for i in 1:nznx; x[i0+i]=inv(rho[i]); end
		elseif(attrib == :rho)
			@inbounds for i in 1:nznx; x[i0+i]=rho[i]; end
		elseif(attrib == :vp)
			@inbounds for i in 1:nznx; x[i0+i]=vp[i]; end
		elseif(attrib == :vs)
			@inbounds for i in 1:nznx; x[i0+i]=vs[i]; end
		else
		end
		i0+=nznx
		end
	end
	χ!(rho, mod.ref.rho,1)
	χ!(vp, mod.ref.vp,1)
	χ!(vs, mod.ref.vs,1)
	return x
end


=#
