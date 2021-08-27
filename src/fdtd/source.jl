
function fill_wavelets!(ipw, iss, wavelets, srcwav, sflags)

	nt = size(wavelets,1)
	sfields=names(srcwav[ipw][iss].d)[1]
	ns=srcwav[ipw][iss][:ns]
	δt = step(srcwav[ipw][iss].grid)

	for sfield in sfields
		for is=1:ns
		snt = length(srcwav[ipw][1].grid)
		if(sflags[ipw] == 0)
			# nothing # just put zeros, no sources added
			for it=1:nt
				wavelets[it][sfield][is] = 0.0
			end
		elseif(sflags[ipw] == 1)
			"ϕ[t] = s[t]"
			for it=1:snt
				source_term = srcwav[ipw][iss].d[sfield][it,is]
				wavelets[it][sfield][is] = source_term
			end
		elseif(sflags[ipw] == -1)
			"source sink for 1"
			"ϕ[t] = -s[t]"
			for it=2:snt
				source_term = srcwav[ipw][iss].d[sfield][it,is]
				wavelets[nt-it+2][sfield][is] = -1.0*source_term
			end
		elseif(sflags[ipw] == 2)
			"ϕ[t] = ∑₁ᵗ⁻¹ s[t]"
			source_term_stack = 0.0;
			if(sfield == :p)
				for it=1:snt-1
					source_term_stack += (srcwav[ipw][iss].d[sfield][it,is] .* δt)
					wavelets[it+1][sfield][is] = source_term_stack
				end
			else
				for it=2:snt-1
					source_term_stack += (((srcwav[ipw][iss].d[sfield][it,is] .* δt) +
					   (srcwav[ipw][iss].d[sfield][it-1,is] .* δt)) * 0.5)
					wavelets[it+1][sfield][is] = source_term_stack
				end

			end
			if(nt > snt)
				wavelets[snt+1:end][sfield][is] = wavelets[snt][sfield][is]
			end
		elseif(sflags[ipw] == -2)
			"source sink for 2"
			# multiplication with -1 for subtraction #time reversal
			# as the source wavelet has to be subtracted before the propagation step, I shift here by one sample"
			source_term_stack = 0.0;
			for it=1:snt-1
				source_term_stack += (srcwav[ipw][iss].d[sfield][it,is] .* δt)
				wavelets[nt-it+1][sfield][is] = -1.0 * source_term_stack
			end
			if(nt > snt)
				nt_diff = nt-snt
				wavelets[1:nt_diff+1][sfield][is] = wavelets[nt_diff+2][sfield][is]
			end
		end
		end
	end
end

struct Source_B1 end
struct Source_B0 end

# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function add_source!(it::Int64, issp::Int64, iss::Int64, pac, pap, ::Source_B1)
	# aliases
	ageom=pac.ageom
	"""
	adding source to pressure field at [it] 
	"""
	for ipw in pac.activepw
		p=pap[ipw].w2[:t];
		wavelets=pap[ipw].ss[issp].wavelets
		isx1=pap[ipw].ss[issp].sindices[:x1]
		isx2=pap[ipw].ss[issp].sindices[:x2]
		isz1=pap[ipw].ss[issp].sindices[:z1]
		isz2=pap[ipw].ss[issp].sindices[:z2]
		ssprayw=pap[ipw].ss[issp].ssprayw
		sfields=pac.sfields[ipw]
		isfields=pac.isfields[ipw]
	if(pac.sflags[ipw] ≠ 0) # only if sflags is non-zero
		for sfield in sfields
		@simd for is = 1:ageom[ipw][iss].ns
			"""
			use wavelets at [it], i.e., sum of source terms
			until [it-1]
			division of source term with δx and δz (see Jan's fdelmodc manual)
			"""
			source_term = wavelets[it][sfield][is] * pac.fc[:δt] * pac.fc[:δxI] * pac.fc[:δzI]
			
			"""
			multiplication with modttI
			"""
			p[sfield][isz1[is], isx1[is]] += 
				source(source_term,ssprayw[1,is], pac, isz1[is], isx1[is],eval(sfield)())
			p[sfield][isz1[is], isx2[is]] += 
				source(source_term,ssprayw[2,is], pac, isz1[is], isx2[is],eval(sfield)())
			p[sfield][isz2[is], isx1[is]] += 
				source(source_term,ssprayw[3,is], pac, isz2[is], isx1[is],eval(sfield)())
			p[sfield][isz2[is], isx2[is]] += 
				source(source_term,ssprayw[4,is], pac, isz2[is], isx2[is],eval(sfield)())
		end
		end
	end
	end
end


# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function add_source!(it::Int64, issp::Int64, iss::Int64, pac, pap, ::Source_B0)
	# aliases
	ageom=pac.ageom
	"""
	adding source to pressure field at [it] 
	"""
	for ipw in pac.activepw
		p=pap[ipw].w2[:t];
		wavelets=pap[ipw].ss[issp].wavelets
		sfields=pac.sfields[ipw]
		isfields=pac.isfields[ipw]
		isx1=pap[ipw].ss[issp].sindices[:x1]
		isz1=pap[ipw].ss[issp].sindices[:z1]
		ssprayw=pap[ipw].ss[issp].ssprayw
	if(pac.sflags[ipw] ≠ 0) # only if sflags is non-zero
		for sfield in sfields
		@simd for is = 1:ageom[ipw].ns[iss]
			"""
			use wavelets at [it], i.e., sum of source terms
			until [it-1]
			division of source term with δx and δz (see Jan's fdelmodc manual)
			"""
			source_term = wavelets[it][sfield][is] * pac.fc[:δt] * pac.fc[:δxI] * pac.fc[:δzI]
			
			"""
			multiplication with modttI
			"""
			p[sfield][isz1[is], isx1[is]] += source(source_term,1.0,pac,isz1[is],isx1[is],eval(sfield)())
		end
	end
	end
	end
end


# on pressure grid
source(source_term, spray, pac, iz, ix, ::p) = source_term * spray * pac.mod[:ttI][iz,ix]

# on vx grid
source(source_term, spray, pac, iz, ix, ::vx) = source_term * spray * pac.mod[:rrvx][iz,ix]

# on vz grid
source(source_term, spray, pac, iz, ix, ::vz) = source_term * spray * pac.mod[:rrvz][iz,ix]



