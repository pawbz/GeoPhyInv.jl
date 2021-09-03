
struct Receiver_B0 end
struct Receiver_B1 end


# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function record!(it::Int64, issp::Int64, iss::Int64, pac, pap, ::Receiver_B1)

	for ipw in pac.activepw
		rinterpolatew=pap[ipw].ss[issp].rinterpolatew
		irx1=pap[ipw].ss[issp].rindices[:x1]
		irx2=pap[ipw].ss[issp].rindices[:x2]
		irz1=pap[ipw].ss[issp].rindices[:z1]
		irz2=pap[ipw].ss[issp].rindices[:z2]
		for rfield in pac.rfields
			recs=pap[ipw].ss[issp].records[rfield]
			pw=pap[ipw].w1[:t][rfield]
			@simd for ir = 1:pac.ageom[ipw][iss].nr
				recs[it,ir]= 
				(
				pw[irz1[ir],irx1[ir]]*
				rinterpolatew[1,ir]+
				pw[irz1[ir],irx2[ir]]*
				rinterpolatew[2,ir]+
				pw[irz2[ir],irx1[ir]]*
				rinterpolatew[3,ir]+
				pw[irz2[ir],irx2[ir]]*
				rinterpolatew[4,ir]
				)
		end
	end
	end
end



# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function record!(it::Int64, issp::Int64, iss::Int64, pac, pap, ::Receiver_B0)
	for ipw in pac.activepw

		rinterpolatew=pap[ipw].ss[issp].rinterpolatew
		irx1=pap[ipw].ss[issp].rindices[:x1]
		irz1=pap[ipw].ss[issp].rindices[:z1]
		for rfield in pac.rfields
			recs=pap[ipw].ss[issp].records[rfield]
			pw=pap[ipw].p.w1[:t][rfield]
			@simd for ir = 1:pac.ageom[ipw][iss].nr
				recs[it,ir]=(pw[irz1[ir],irx1[ir]])
		end
	end
	end
end


