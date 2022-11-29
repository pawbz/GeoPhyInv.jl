
# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function record!(it::Int64, issp::Int64, iss::Int64, pac, pap)

    for ipw in pac.activepw
        rinterpolatew = pap[ipw].ss[issp].rinterpolatew
        for rfield in pac.rfields
            pw = view(pap[ipw].w1[:t][rfield],:)
            recs = pap[ipw].ss[issp].records[rfield][it]
            w=rinterpolatew[rfield]
            mul!(recs, w, pw)
        end
    end
end


