"""
* `g` : input gradients w.r.t. KI and rhoI
* `gout` : output gradients according to parameterization
"""
function chainrule!(gout::AbstractVector, g::AbstractVector, mod::Medium, parameterization)
	nznx=prod(length.(mod.mgrid))
	(length(gout)≠(count(parameterization.≠ :null)*nznx)) &&  error("size x")

	rho=mod[:rho]; 
	vp=mod[:vp]; 
	# vs=mod[:vs]; # not implemented

	if(parameterization == [:χKI, :χrhoI, :null]) 
		@inbounds for i in 1:nznx
			gout[i]=χg(g[i],mod.ref[:KI],1)
			gout[nznx+i]=χg(g[nznx+i],mod.ref[:rhoI],1)
		end
	elseif(parameterization == [:χKI, :null, :null]) 
		@inbounds for i in 1:nznx
			gout[i]=χg(g[i],mod.ref[:KI],1)
		end
	elseif(parameterization == [:χvp, :χrho, :null]) 
		@inbounds for i in 1:nznx
			# derivative w.r.t vp
			gout[i]= χg(grad_of_vp(g[i], vp[i], rho[i]), mod.ref[:vp],1)
			# derivative w.r.t. rho
			gout[nznx+i]=χg(grad_of_rho(g[i],g[nznx+i],vp[i],rho[i]), mod.ref[:rho],1)
		end
	elseif(parameterization == [:χvp, :null, :null]) 
		@inbounds for i in 1:nznx
			gout[i]= χg(grad_of_vp(g[i], vp[i], rho[i]), mod.ref[:vp],1)
		end
	elseif(parameterization == [:null, :χrho, :null]) 
		@inbounds for i in 1:nznx
			gout[i]=χg(grad_of_rho(g[i],vp[i],rho[i]), mod.ref[:rho],1)
		end
	end

end



grad_of_vp(gKI, vp, rho) = -2.0 * gKI * inv(vp*vp*vp*rho)
grad_of_rho(gKI, grhoI, vp, rho) = -1.0 * inv(abs2(rho)) * (inv(vp*vp)*gKI + grhoI) 
grad_of_rho(grhoI, vp, rho) = -1.0 * inv(abs2(rho)) * (grhoI) 

"""
No different from the previous case, but...
"""
function pert_chainrule!(gout::AbstractVector, g::AbstractVector, mod::Medium, parameterization)
	nznx=prod(length.(mod.mgrid))
	(length(gout)≠(count(parameterization.≠ :null)*nznx)) &&  error("size x")
	fill!(gout, 0.0)
	if(parameterization == [:χKI, :χrhoI, :null]) 
		@inbounds for i in 1:nznx
			gout[i]=g[i]*mod.ref[:KI]
			gout[nznx+i]=g[nznx+i]*mod.ref[:rhoI]
		end
	elseif(parameterization == [:χKI, :null, :null]) 
		@inbounds for i in 1:nznx
			gout[i]=χg(g[i],mod.ref[:KI],1)
		end
	elseif(parameterization == [:χvp, :χrho, :null]) 
		@inbounds for i in 1:nznx
			# derivative w.r.t vp
			gout[i]= -2.0*mod.ref[:KI]*g[i]
			# derivative w.r.t. rho
			gout[nznx+i]=-1.0*mod.ref[:KI]*g[i]-1.0*mod.ref[:rhoI]*g[nznx+i]
		end
	elseif(parameterization == [:χvp, :null, :null]) 
		@inbounds for i in 1:nznx
			gout[i]=-2.0*mod.ref[:KI]*g[i]
		end
	elseif(parameterization == [:null, :χrho, :null]) 
		@inbounds for i in 1:nznx
			gout[i]=-1.0*mod.ref[:rhoI]*g[i]
		end
	end
end









"""
Use chain rule to output gradients with 
respect to χvp and χrho from  gradients 
with respect to KI and rhoI.

# Arguments

* `gmod::Seismic` : gradient model
* `mod::Seismic` : model required for chain rule
* `g1` : gradient of an objective function with respect `attribs[1]`
* `g2` : gradient of an objective function with respect `attribs[2]`
* `attribs::Vector{Symbol}=[:χKI, :χrhoI]` :  
* `flag::Int64=1` :
  * `=1` updates `gmod` using `g1` and `g2`
  * `=-1` updates `g1` and `g2` using `gmod`
"""
function chainrule!(
		      gmod::Medium,
		      mod::Medium,
		      g,
		      attribvec::Vector{Symbol}=[:χKI, :χrhoI, :null],
		      flag::Int64=1
		      )

	nznx=prod(length.(mod.mgrid))
	(length(g)≠(count(attribvec.≠ :null)*nznx)) &&  error("size x")

	rho=mod[:rho]; vp=mod[:vp]; 
	# vs=mod[:vs];
	if(flag == 1)
		if(attribvec == [:χKI, :χrhoI, :null]) 
			@inbounds for i in 1:nznx
				g[i]=χg(g[i],mod.ref[:KI],-1)
				g[nznx+i]=χg(g[nznx+i],mod.ref[:rhoI],-1)
				# gradient w.r.t  χrho
				gmod.m[:χrho][i]= -1.0*abs2(inv(rho[i]))*g[nznx+i]-inv(vp[i]*vp[i]*rho[i]*rho[i])*g[i]
				# gradient w.r.t χvp
				gmod.m[:χvp][i] = -2.0*inv(vp[i]*vp[i]*vp[i])*inv(rho[i])*g[i]
				g[i]=χg(g[i],mod.ref[:KI],1)
				g[nznx+i]=χg(g[nznx+i],mod.ref[:rhoI],1)
			end
			χg!(gmod.m[:χrho],mod.ref[:rho],1)
			χg!(gmod.m[:χvp],mod.ref[:vp],1)
		elseif(attribvec == [:χKI, :null, :null]) 
			@inbounds for i in 1:nznx
				g[i]=χg(g[i],mod.ref[:KI],-1)
				# gradient w.r.t χvp
				gmod.m[:χvp][i] = -2.0*inv(vp[i]*vp[i]*vp[i])*inv(rho[i])*g[i]
				g[i]=χg(g[i],mod.ref[:KI],1)
			end
			χg!(gmod.m[:χvp],mod.ref[:vp],1)

		elseif(attribvec == [:χvp, :χrho, :null]) 
			@inbounds for i in 1:nznx; gmod.m[:χvp][i]=g[i]; end
			@inbounds for i in 1:nznx; gmod.m[:χrho][i]=g[nznx+i]; end
		elseif(attribvec == [:null, :χrho, :null]) 
			@inbounds for i in 1:nznx; gmod.m[:χrho][i]=g[i]; end
		elseif(attribvec == [:χvp, :null, :null]) 
			@inbounds for i in 1:nznx; gmod.m[:χvp][i]=g[i]; end
		else
			error("invalid attribs")
		end
	elseif(flag == -1)
		if(attribvec == [:χKI, :χrhoI, :null])
			χg!(gmod.m[:χvp], mod.ref[:vp], -1) # undo it later
			χg!(gmod.m[:χrho], mod.ref[:rho], -1) # undo it later
			gvp = gmod.m[:χvp]
			grho = gmod.m[:χrho]
			@inbounds for i in 1:nznx
				g[i] = -0.5 * mod.ref[:KI] * gvp[i] * rho[i] *  vp[i]^3
			end
			@inbounds for i in 1:nznx
				g[nznx+i] = (grho[i] * mod.ref[:rhoI] * -1. * rho[i]*rho[i]) + 
				(0.5 * mod.ref[:rhoI] .* gvp[i] .*  vp[i] .* rho[i])
			end
			χg!(gvp,mod.ref[:vp],1)
			χg!(grho,mod.ref[:rho],1)
		elseif(attribvec == [:χKI, :null, :null])
			χg!(gmod.m[:χvp], mod.ref[:vp], -1) # undo it later
			gvp = gmod.m[:χvp]
			@inbounds for i in 1:nznx
				g[i] = -0.5 * mod.ref[:KI] * gvp[i] * rho[i] *  vp[i]^3
			end
			χg!(gvp,mod.ref[:vp],1)
		elseif(attribvec == [:χvp, :χrho, :null]) 
			@inbounds for i in 1:nznx
				g[i] = gmod.m[:χvp][i]
				g[nznx+i] = gmod.m[:χrho][i]
			end
		elseif(attribvec == [:χvp, :null, :null]) 
			@inbounds for i in 1:nznx
				g[i] = gmod.m[:χvp][i]
			end
		elseif(attribvec == [:null, :χrho, :null]) 
			@inbounds for i in 1:nznx
				g[i] = gmod.m[:χrho][i]
			end
		else
			error("invalid attribs")
		end
	else
		error("invalid flag")
	end
end
