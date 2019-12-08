
function fill_wavelets!(iss::Int64, wavelets::Array{Array{Float64,2},2}, srcwav::Array{SrcWav}, sflags::Vector{Int64})

	npw = size(wavelets,1)
	nt = size(wavelets,2)
	δt = step(srcwav[1][1].grid)
	for ipw=1:npw
		ns, snfield = size(wavelets[ipw,1]) # ns may vary with ipw
		for ifield=1:snfield, is=1:ns
			snt = length(srcwav[ipw][1].grid)
			if(sflags[ipw] == 0)
				# nothing # just put zeros, no sources added
				for it=1:nt
					wavelets[ipw,it][is,ifield] = 0.0
				end
			elseif(sflags[ipw] == 1)
				"ϕ[t] = s[t]"
				for it=1:snt
					source_term = srcwav[ipw][iss].d[ifield][it,is]
					wavelets[ipw,it][is,ifield] = source_term
				end
			elseif(sflags[ipw] == -1)
				"source sink for 1"
				"ϕ[t] = -s[t]"
				for it=2:snt
					source_term = srcwav[ipw][iss].d[ifield][it,is]
					wavelets[ipw,nt-it+2][is,ifield] = -1.0*source_term
				end
			elseif(sflags[ipw] == 2)
				"ϕ[t] = ∑₁ᵗ⁻¹ s[t]"
				source_term_stack = 0.0;
				if(ifield == 1)
					for it=1:snt-1
						source_term_stack += (srcwav[ipw][iss].d[ifield][it,is] .* δt)
						wavelets[ipw,it+1][is,ifield] = source_term_stack
					end
				else
					for it=2:snt-1
						source_term_stack += (((srcwav[ipw][iss].d[ifield][it,is] .* δt) +
						   (srcwav[ipw][iss].d[ifield][it-1,is] .* δt)) * 0.5)
						wavelets[ipw,it+1][is,ifield] = source_term_stack
					end

				end
				if(nt > snt)
					wavelets[ipw,snt+1:end][is,ifield] = wavelets[ipw,snt][is,ifield]
				end
			elseif(sflags[ipw] == -2)
				"source sink for 2"
				# multiplication with -1 for subtraction #time reversal
				# as the source wavelet has to be subtracted before the propagation step, I shift here by one sample"
				source_term_stack = 0.0;
				for it=1:snt-1
					source_term_stack += (srcwav[ipw][iss].d[ifield][it,is] .* δt)
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
@inbounds @fastmath function add_source!(it::Int64, issp::Int64, iss::Int64, pac::PFdtdc, pass::Vector{PFdtdss}, pap::PFdtdp, ::Source_B1)
	# aliases
	p=pap.p;
	wavelets=pass[issp].wavelets
	geom=pac.geom
	isx1=pass[issp].isx1
	isx2=pass[issp].isx2
	isz1=pass[issp].isz1
	isz2=pass[issp].isz2
	ssprayw=pass[issp].ssprayw
	modttI=pac.modttI
	modrr=pac.modrrvz
	"""
	adding source to pressure field at [it] 
	"""
	for ipw in pac.activepw
		sfields=pac.sfields[ipw]
		isfields=pac.isfields[ipw]
	if(pac.sflags[ipw] ≠ 0) # only if sflags is non-zero
		pw=p[ipw]	
		for (iff, ifield) in enumerate(isfields)
		@simd for is = 1:geom[ipw][iss].ns
			"""
			use wavelets at [it], i.e., sum of source terms
			until [it-1]
			division of source term with δx and δz (see Jan's fdelmodc manual)
			"""
			source_term = wavelets[ipw,it][is, iff] * pac.δt * pac.δxI * pac.δzI
			
			"""
			multiplication with modttI
			"""
			pw[isz1[ipw][is], isx1[ipw][is],ifield] += 
				source(source_term,ssprayw[ipw][1,is], pac, isz1[ipw][is], isx1[ipw][is],eval(sfields[iff])())
			pw[isz1[ipw][is], isx2[ipw][is],ifield] += 
				source(source_term,ssprayw[ipw][2,is], pac, isz1[ipw][is], isx2[ipw][is],eval(sfields[iff])())
			pw[isz2[ipw][is], isx1[ipw][is],ifield] += 
				source(source_term,ssprayw[ipw][3,is], pac, isz2[ipw][is], isx1[ipw][is],eval(sfields[iff])())
			pw[isz2[ipw][is], isx2[ipw][is],ifield] += 
				source(source_term,ssprayw[ipw][4,is], pac, isz2[ipw][is], isx2[ipw][is],eval(sfields[iff])())
		end
		end
	end
	end
end


# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function add_source!(it::Int64, issp::Int64, iss::Int64, pac::PFdtdc, pass::Vector{PFdtdss}, pap::PFdtdp, ::Source_B0)
	# aliases
	p=pap.p;
	wavelets=pass[issp].wavelets
	geom=pac.geom
	isx1=pass[issp].isx1
	isz1=pass[issp].isz1
	ssprayw=pass[issp].ssprayw
	modttI=pac.modttI
	"""
	adding source to pressure field at [it] 
	"""
	for ipw in pac.activepw
		sfields=pac.sfields[ipw]
		isfields=pac.isfields[ipw]
	if(pac.sflags[ipw] ≠ 0) # only if sflags is non-zero
		pw=p[ipw]
		for (iff, ifield) in enumerate(isfields)
		@simd for is = 1:geom[ipw].ns[iss]
			"""
			use wavelets at [it], i.e., sum of source terms
			until [it-1]
			division of source term with δx and δz (see Jan's fdelmodc manual)
			"""
			source_term = wavelets[ipw,it][is, iff] * pac.δt * pac.δxI * pac.δzI
			
			"""
			multiplication with modttI
			"""
			pw[isz1[ipw][is], isx1[ipw][is],ifield] += source(source_term,1.0,pac,isz1[ipw][is],isx1[ipw][is],eval(sfields[iff])())
		end
	end
	end
	end
end


# on pressure grid
source(source_term, spray, pac, iz, ix, ::P) = source_term * spray * pac.modttI[iz,ix]

# on Vx grid
source(source_term, spray, pac, iz, ix, ::Vx) = source_term * spray * pac.modrrvx[iz,ix]

# on Vz grid
source(source_term, spray, pac, iz, ix, ::Vz) = source_term * spray * pac.modrrvz[iz,ix]



