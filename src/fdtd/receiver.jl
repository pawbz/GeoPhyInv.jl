
# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function record!(it::Int64, issp::Int64, iss::Int64, pac, pap, activepw, rfields=pac.rfields)
    for ipw in activepw
        rinterpolatew = pap[ipw].ss[issp].rinterpolatew
        rfields1 = intersect(rfields, pac.rfields)
        for rfield in rfields1
            pw = view(pap[ipw].w1[:t][rfield],:)
            recs = pap[ipw].ss[issp].records[rfield][it]
            w=rinterpolatew[rfield]
            mul!(recs, w, pw)
        end
    end
end


function update_datamat!(
    rfield,
    ipw,
    pac::P_common,
    pap::Vector{P_x_worker_x_pw{N,B}},
) where {N,B}
    datamat = pac.datamat
    pass = pap[ipw].ss
    for issp = 1:length(pass)
        iss = pass[issp].iss
        records = pass[issp].records
        for it = 1:pac.ic[:nt]
            r = Array(records[rfield][it])
            d = view(datamat, it, :, iss)
            copyto!(d, r)
        end
    end
end

function update_data!(rfield, ipw, pac::P_common)
    datamat = pac.datamat
    for iss = 1:length(pac.ageom[1])
        data = pac.data[ipw][iss].d[rfield]
        for ir = 1:pac.ageom[ipw][iss].nr
            for it = 1:pac.ic[:nt]
                data[it, ir] = datamat[it, ir, iss]
            end
        end
    end
end

