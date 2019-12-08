
struct Receiver_B0 end
struct Receiver_B1 end


# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function record!(it::Int64, issp::Int64, iss::Int64, pac::Fdtdc, pass::Vector{Fdtdss}, pap::Fdtdp, ::Receiver_B1)
	p=pap.p
	rinterpolatew=pass[issp].rinterpolatew
	irx1=pass[issp].irx1
	irx2=pass[issp].irx2
	irz1=pass[issp].irz1
	irz2=pass[issp].irz2

	for ipw in pac.activepw
		pw=p[ipw]
		recs=pass[issp].records[ipw]
		for (ifieldr, ifield) in enumerate(pac.irfields)
			@simd for ir = 1:pac.acqgeom[ipw][iss].nr
				recs[it,ir,ifieldr]= 
				(
				pw[irz1[ipw][ir],irx1[ipw][ir],ifield]*
				rinterpolatew[ipw][1,ir]+
				pw[irz1[ipw][ir],irx2[ipw][ir],ifield]*
				rinterpolatew[ipw][2,ir]+
				pw[irz2[ipw][ir],irx1[ipw][ir],ifield]*
				rinterpolatew[ipw][3,ir]+
				pw[irz2[ipw][ir],irx2[ipw][ir],ifield]*
				rinterpolatew[ipw][4,ir]
				)
		end
	end
	end
end



# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function record!(it::Int64, issp::Int64, iss::Int64, pac::Fdtdc, pass::Vector{Fdtdss}, pap::Fdtdp, ::Receiver_B0)
	p=pap.p
	rinterpolatew=pass[issp].rinterpolatew
	irx1=pass[issp].irx1
	irz1=pass[issp].irz1

	for ipw in pac.activepw
		pw=p[ipw]
		recs=pass[issp].records[ipw]
		for (ifieldr, ifield) in enumerate(pac.irfields)
			@simd for ir = 1:pac.acqgeom[ipw][iss].nr
				recs[it,ir,ifieldr]=(pw[irz1[ipw][ir],irx1[ipw][ir],ifield])
		end
	end
	end
end


