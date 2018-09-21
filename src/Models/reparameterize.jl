
"""
Re-parameterization routine 
that modifies the fields 
`χvp` and `χρ` of an input seismic model
using two input vectors.

# Arguments

* `mod::Seismic` : to be updated
* `x1::Array{Float64,2}` : contrast of inverse bulk modulus
* `x2::Array{Float64,2}` : contrast of inverse density
* `attribvec:::Vector{Symbol}` : [:χKI, :χρI]
"""
function Seismic_reparameterize!(
	              mod::Seismic,
		      x,
		      attribvec::Vector{Symbol}=[:χKI, :χρI, :null]
		      )
	iszero(mod) ? error("mod cannot be zero") : nothing
	nznx=prod(length.(mod.mgrid))
	(length(x)≠(count(attribvec.≠ :null)*nznx)) &&  error("size x")
	if(attribvec == [:χKI, :χρI, :null]) 
		@inbounds for i in 1:nznx 
			K=inv(χ(x[i],mod.ref.KI,-1))
			ρ=inv(χ(x[nznx+i],mod.ref.ρI,-1))
			mod.χvp[i]=χ(sqrt(K*inv(ρ)),mod.ref.vp,1)
			mod.χρ[i]=χ(ρ,mod.ref.ρ,1)
		end
	elseif(attribvec == [:χKI, :null, :null]) 
		@inbounds for i in 1:nznx 
			K=inv(χ(x[i],mod.ref.KI,-1))
			ρ=χ(mod.χρ[i],mod.ref.ρ,-1)
			mod.χvp[i]=χ(sqrt(K*inv(ρ)),mod.ref.vp,1)
			mod.χρ[i]=χ(ρ,mod.ref.ρ,1)
		end
	elseif(attribvec == [:χvp, :χρ, :null]) 
		@inbounds for i in 1:nznx; mod.χvp[i]=x[i]; end
		@inbounds for i in 1:nznx; mod.χρ[i]=x[nznx+i]; end
	elseif(attribvec == [:null, :χρ, :null]) 
		@inbounds for i in 1:nznx; mod.χρ[i]=x[i]; end
	elseif(attribvec == [:χvp, :null, :null]) 
		@inbounds for i in 1:nznx; mod.χvp[i]=x[i]; end
	else
		error("invalid attribvec")
	end
	return mod
end




"""
Input perturbations corresponding to `parameterization`.
Output perturbations of KI and rhoI, i.e., without the contrast.
`mod` is just used for reference values.
"""
function pert_reparameterize!(δxout, δx, mod, parameterization)
	nznx=prod(length.(mod.mgrid))
	fill!(δxout, 0.0)
	if(parameterization == [:χKI, :χρI, :null]) 
		@inbounds for i in 1:nznx
			δxout[i]=δx[i]*mod.ref.KI
			δxout[nznx+i]=δx[nznx+i]*mod.ref.ρI
		end
	elseif(parameterization == [:χKI, :null, :null]) 
		@inbounds for i in 1:nznx
			δxout[i]=δx[i]*mod.ref.KI
		end
	elseif(parameterization == [:χvp, :χρ, :null]) 
		@inbounds for i in 1:nznx
			δxout[i]=(-1.0*δx[nznx+i]-2.0*δx[i])*mod.ref.KI
			δxout[nznx+i]=(-1.0*δx[nznx+i])*mod.ref.ρI
		end
	elseif(parameterization == [:χvp, :null, :null]) 
		@inbounds for i in 1:nznx
			δxout[i]=-2.0*δx[i]*mod.ref.KI
		end
	elseif(parameterization == [:null, :χρ, :null]) 
		@inbounds for i in 1:nznx
			δxout[i]=-1.0*δx[i]*mod.ref.ρI
		end
	end
end


#pert_of_KI(dvp, dρ, vp0, ρ0) = - 2. / (ρ0 * vp0^3) * dvp - 1 / (ρ0^2 * vp0^2) * dρ
#pert_of_ρI(dρ, vp0, ρ0) = -1. / (ρ0^2) * dρ
