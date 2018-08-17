
# the input derivatives are w.r.t KI=

"""
* `g` : input gradients w.r.t. KI and ρI
* `gout` : output gradients according to parameterization
"""
function gradient_chainrule!(gout, g, mod, parameterization)
	nznx=mod.mgrid.nz*mod.mgrid.nx
	(length(gout)≠(count(parameterization.≠ :null)*nznx)) &&  error("size x")

	ρ=mod.χρ; χ!(ρ, mod.ref.ρ,-1) # undo it later
	vp=mod.χvp; χ!(vp, mod.ref.vp,-1) # undo it later
	vs=mod.χvs; χ!(vs, mod.ref.vs,-1) # undo it later

	if(parameterization == [:χKI, :χρI, :null]) 
		@inbounds for i in 1:nznx
			gout[i]=χg(g[i],mod.ref.KI,1)
			gout[nznx+i]=χg(g[nznx+i],mod.ref.ρI,1)
		end
	elseif(parameterization == [:χKI, :null, :null]) 
		@inbounds for i in 1:nznx
			gout[i]=χg(g[i],mod.ref.KI,1)
		end
	elseif(parameterization == [:χvp, :χρ, :null]) 
		@inbounds for i in 1:nznx
			# derivative w.r.t vp
			gout[i]= χg(grad_of_vp(g[i], vp[i], ρ[i]), mod.ref.vp,1)
			# derivative w.r.t. ρ
			gout[nznx+i]=χg(grad_of_ρ(g[i],g[nznx+i],vp[i],ρ[i]), mod.ref.ρ,1)
		end
	elseif(parameterization == [:χvp, :null, :null]) 
		@inbounds for i in 1:nznx
			gout[i]= χg(grad_of_vp(g[i], vp[i], ρ[i]), mod.ref.vp,1)
		end
	elseif(parameterization == [:null, :χρ, :null]) 
		@inbounds for i in 1:nznx
			gout[i]=χg(grad_of_ρ(g[i],g[nznx+i],vp[i],ρ[i]), mod.ref.ρ,1)
		end
	end

	χ!(ρ, mod.ref.ρ,1) 
	χ!(vp, mod.ref.vp,1)
	χ!(vs, mod.ref.vs,1)

end



grad_of_vp(gKI, vp, ρ) = -2. * gKI * inv(vp*vp*vp*ρ)
grad_of_ρ(gKI, gρI, vp, ρ) = -1. * inv(abs2(ρ)) * (inv(vp*vp)*gKI + gρI) 








"""
Use chain rule to output gradients with 
respect to χvp and χρ from  gradients 
with respect to KI and ρI.

# Arguments

* `gmod::Seismic` : gradient model
* `mod::Seismic` : model required for chain rule
* `g1` : gradient of an objective function with respect `attribs[1]`
* `g2` : gradient of an objective function with respect `attribs[2]`
* `attribs::Vector{Symbol}=[:χKI, :χρI]` :  
* `flag::Int64=1` :
  * `=1` updates `gmod` using `g1` and `g2`
  * `=-1` updates `g1` and `g2` using `gmod`
"""
function Seismic_chainrule!(
		      gmod::Seismic,
		      mod::Seismic,
		      g,
		      attribvec::Vector{Symbol}=[:χKI, :χρI, :null],
		      flag::Int64=1
		      )

	nznx=mod.mgrid.nz*mod.mgrid.nx
	(length(g)≠(count(attribvec.≠ :null)*nznx)) &&  error("size x")

	ρ=mod.χρ; χ!(ρ, mod.ref.ρ,-1) # undo it later
	vp=mod.χvp; χ!(vp, mod.ref.vp,-1) # undo it later
	vs=mod.χvs; χ!(vs, mod.ref.vs,-1) # undo it later
	if(flag == 1)
		if(attribvec == [:χKI, :χρI, :null]) 
			@inbounds for i in 1:nznx
				g[i]=χg(g[i],mod.ref.KI,-1)
				g[nznx+i]=χg(g[nznx+i],mod.ref.ρI,-1)
				# gradient w.r.t  χρ
				gmod.χρ[i]= -1.0*abs2(inv(ρ[i]))*g[nznx+i]-inv(vp[i]*vp[i]*ρ[i]*ρ[i])*g[i]
				# gradient w.r.t χvp
				gmod.χvp[i] = -2.0*inv(vp[i]*vp[i]*vp[i])*inv(ρ[i])*g[i]
				g[i]=χg(g[i],mod.ref.KI,1)
				g[nznx+i]=χg(g[nznx+i],mod.ref.ρI,1)
			end
			χg!(gmod.χρ,mod.ref.ρ,1)
			χg!(gmod.χvp,mod.ref.vp,1)
		elseif(attribvec == [:χKI, :null, :null]) 
			@inbounds for i in 1:nznx
				g[i]=χg(g[i],mod.ref.KI,-1)
				# gradient w.r.t χvp
				gmod.χvp[i] = -2.0*inv(vp[i]*vp[i]*vp[i])*inv(ρ[i])*g[i]
				g[i]=χg(g[i],mod.ref.KI,1)
			end
			χg!(gmod.χvp,mod.ref.vp,1)

		elseif(attribvec == [:χvp, :χρ, :null]) 
			@inbounds for i in 1:nznx; gmod.χvp[i]=g[i]; end
			@inbounds for i in 1:nznx; gmod.χρ[i]=g[nznx+i]; end
		elseif(attribvec == [:null, :χρ, :null]) 
			@inbounds for i in 1:nznx; gmod.χρ[i]=g[i]; end
		elseif(attribvec == [:χvp, :null, :null]) 
			@inbounds for i in 1:nznx; gmod.χvp[i]=g[i]; end
		else
			error("invalid attribs")
		end
	elseif(flag == -1)
		if(attribvec == [:χKI, :χρI, :null])
			χg!(gmod.χvp, mod.ref.vp, -1) # undo it later
			χg!(gmod.χρ, mod.ref.ρ, -1) # undo it later
			gvp = gmod.χvp
			gρ = gmod.χρ
			@inbounds for i in 1:nznx
				g[i] = -0.5 * mod.ref.KI * gvp[i] * ρ[i] *  vp[i]^3
			end
			@inbounds for i in 1:nznx
				g[nznx+i] = (gρ[i] * mod.ref.ρI * -1. * ρ[i]*ρ[i]) + 
					(0.5 * mod.ref.ρI .* gvp[i] .*  vp[i] .* ρ[i])
			end
			χg!(gvp,mod.ref.vp,1)
			χg!(gρ,mod.ref.ρ,1)
		elseif(attribvec == [:χKI, :null, :null])
			χg!(gmod.χvp, mod.ref.vp, -1) # undo it later
			gvp = gmod.χvp
			@inbounds for i in 1:nznx
				g[i] = -0.5 * mod.ref.KI * gvp[i] * ρ[i] *  vp[i]^3
			end
			χg!(gvp,mod.ref.vp,1)
		elseif(attribvec == [:χvp, :χρ, :null]) 
			@inbounds for i in 1:nznx
				g[i] = gmod.χvp[i]
				g[nznx+i] = gmod.χρ[i]
			end
		elseif(attribvec == [:χvp, :null, :null]) 
			@inbounds for i in 1:nznx
				g[i] = gmod.χvp[i]
			end
		elseif(attribvec == [:null, :χρ, :null]) 
			@inbounds for i in 1:nznx
				g[i] = gmod.χρ[i]
			end
		else
			error("invalid attribs")
		end
	else
		error("invalid flag")
	end
	χ!(ρ, mod.ref.ρ,1) 
	χ!(vp, mod.ref.vp,1)
	χ!(vs, mod.ref.vs,1)
end
