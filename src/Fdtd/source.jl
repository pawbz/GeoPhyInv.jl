
function fill_wavelets!(iss::Int64, wavelets::Array{Array{Float64,2},2}, acqsrc::Array{Acquisition.Src}, sflags::Vector{Int64})

	npw = size(wavelets,1)
	nt = size(wavelets,2)
	δt = acqsrc[1].tgrid.δx
	for ipw=1:npw
		ns, snfield = size(wavelets[ipw,1]) # ns may vary with ipw
		for ifield=1:snfield, is=1:ns
			snt = acqsrc[ipw].tgrid.nx;
			if(sflags[ipw] == 0)
				# nothing # just put zeros, no sources added
				for it=1:nt
					wavelets[ipw,it][is,ifield] = 0.0
				end
			elseif(sflags[ipw] == 1)
				"ϕ[t] = s[t]"
				for it=1:snt
					source_term = acqsrc[ipw].wav[iss,ifield][it,is]
					wavelets[ipw,it][is,ifield] = source_term
				end
			elseif(sflags[ipw] == -1)
				"ϕ[t] = s[t]"
				for it=1:snt
					source_term = acqsrc[ipw].wav[iss,ifield][it,is]
					wavelets[ipw,nt-it+1][is,ifield] = -1.0*source_term
				end
			elseif(sflags[ipw] == 2)
				"ϕ[t] = ∑₁ᵗ⁻¹ s[t]"
				source_term_stack = 0.0;
				if(ifield == 1)
					for it=1:snt-1
						source_term_stack += (acqsrc[ipw].wav[iss,ifield][it,is] .* δt)
						wavelets[ipw,it+1][is,ifield] = source_term_stack
					end
				else
					for it=2:snt-1
						source_term_stack += (((acqsrc[ipw].wav[iss,ifield][it,is] .* δt) +
						   (acqsrc[ipw].wav[iss,ifield][it-1,is] .* δt)) * 0.5)
						wavelets[ipw,it+1][is,ifield] = source_term_stack
					end

				end
				if(nt > snt)
					wavelets[ipw,snt+1:end][is,ifield] = wavelets[ipw,snt][is,ifield]
				end
			elseif(sflags[ipw] == 3)
				"use this to add source sink: need during adjoint propagation from boundary"
				"multiplication with -1 for subtraction"
				"time reversal"
				"as the source wavelet has to be subtracted before the propagation step, I shift here by one sample"
				"ϕ[t] = "
				source_term_stack = 0.0;
				for it=1:snt-1
					source_term_stack += (acqsrc[ipw].wav[iss,ifield][it,is] .* δt)
					wavelets[ipw,nt-it+1][is,ifield] = -1.0 * source_term_stack
				end
				if(nt > snt)
					nt_diff = nt-snt
					wavelets[ipw,1:nt_diff+1][is,ifield] = wavelets[ipw,nt_diff+2][is,ifield]
				end
			end
		end
	end
end

struct Source_B1 end
struct Source_B0 end

# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function add_source!(it::Int64, issp::Int64, iss::Int64, pac::Paramc, pass::Vector{Paramss}, pap::Paramp, ::Source_B1)
	# aliases
	p=pap.p;
	wavelets=pass[issp].wavelets
	acqgeom=pac.acqgeom
	isx1=pass[issp].isx1
	isx2=pass[issp].isx2
	isz1=pass[issp].isz1
	isz2=pass[issp].isz2
	ssprayw=pass[issp].ssprayw
	modttI=pac.modttI
	"""
	adding source to pressure field at [it] 
	"""
	for ipw in pac.activepw
	if(pac.sflags[ipw] ≠ 0) # only if sflags is non-zero
		pw=p[ipw]	
		for (ifields, ifield) in enumerate(pac.isfields[ipw])
		@simd for is = 1:acqgeom[ipw].ns[iss]
			"""
			use wavelets at [it], i.e., sum of source terms
			until [it-1]
			division of source term with δx and δz (see Jan's fdelmodc manual)
			"""
			source_term = wavelets[ipw,it][is, ifields] * pac.δt * pac.δxI * pac.δzI

						if(it==1)
							println("source term t=0*^^***** ", source_term)
						end
			
			"""
			multiplication with modttI
			"""
			pw[isz1[ipw][is], isx1[ipw][is],ifield] += 
				source_term * 
				ssprayw[ipw][1,is] * 
				modttI[isz1[ipw][is], isx1[ipw][is]]  
			pw[isz1[ipw][is], isx2[ipw][is],ifield] += 
				source_term * 
				ssprayw[ipw][2,is] * 
				modttI[isz1[ipw][is], isx2[ipw][is]]
			pw[isz2[ipw][is], isx1[ipw][is],ifield] += 
				source_term * 
				ssprayw[ipw][3,is] * 
				modttI[isz2[ipw][is], isx1[ipw][is]]
			pw[isz2[ipw][is], isx2[ipw][is],ifield] += 
				source_term * 
				ssprayw[ipw][4,is] * 
				modttI[isz2[ipw][is], isx2[ipw][is]]
		end
		end
	end
	end
end



# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function add_source!(it::Int64, issp::Int64, iss::Int64, pac::Paramc, pass::Vector{Paramss}, pap::Paramp, ::Source_B0)
	# aliases
	p=pap.p;
	wavelets=pass[issp].wavelets
	acqgeom=pac.acqgeom
	isx1=pass[issp].isx1
	isz1=pass[issp].isz1
	ssprayw=pass[issp].ssprayw
	modttI=pac.modttI
	"""
	adding source to pressure field at [it] 
	"""
	for ipw in pac.activepw
	if(pac.sflags[ipw] ≠ 0) # only if sflags is non-zero
		pw=p[ipw]
		for (ifields, ifield) in enumerate(pac.isfields[ipw])
		@simd for is = 1:acqgeom[ipw].ns[iss]
			"""
			use wavelets at [it], i.e., sum of source terms
			until [it-1]
			division of source term with δx and δz (see Jan's fdelmodc manual)
			"""
			source_term = wavelets[ipw,it][is, ifields] * pac.δt * pac.δxI * pac.δzI
			
			"""
			multiplication with modttI
			"""
			pw[isz1[ipw][is], isx1[ipw][is],ifield] += source_term * modttI[isz1[ipw][is], isx1[ipw][is]]  
		end
	end
	end
	end
end


