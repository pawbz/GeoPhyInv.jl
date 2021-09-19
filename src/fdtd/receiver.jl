
struct Receiver_B0 end
struct Receiver_B1 end


# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function record!(
    it::Int64,
    issp::Int64,
    iss::Int64,
    pac,
    pap,
)

    for ipw in pac.activepw
        rinterpolatew = pap[ipw].ss[issp].rinterpolatew
        rindices = pap[ipw].ss[issp].rindices
        for rfield in pac.rfields
            recs = pap[ipw].ss[issp].records[rfield]
            pw = pap[ipw].w1[:t][rfield]

            @simd for ir = 1:pac.ageom[ipw][iss].nr
                recs[it, ir] = 0.0
                for (i, ri) in enumerate(rindices[ir])
                    recs[it, ir] += pw[ri] * rinterpolatew[ir][i]
                end
            end
        end
    end
end



