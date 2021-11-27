"""
Use input srcwav to fill wavelets and returns frequency bounds.
"""
function fill_wavelets!(ipw, iss, wavelets, srcwav, sflags)

    nt = size(wavelets, 1)
    sfields = names(srcwav[ipw][iss].d)[1]
    ns = srcwav[ipw][iss][:ns]
    δt = step(srcwav[ipw][iss].grid)

    freqmin = 0.0
    freqpeaks = []
    freqmax = Inf
    for sfield in sfields
        snt = length(srcwav[ipw][1].grid)
        if (sflags[ipw] == 0)
            # nothing # just put zeros, no sources added
            for it = 1:nt
                w = zeros(ns)
                wavelets[it][sfield] = Data.Array(w)
            end
        elseif (sflags[ipw] == 1)
            "ϕ[t] = s[t]"
            for it = 1:snt
                w = [srcwav[ipw][iss].d[sfield][it, is] for is = 1:ns]
                wavelets[it][sfield] = Data.Array(w)
            end
        elseif (sflags[ipw] == -1)
            "source sink for 1"
            "ϕ[t] = -s[t]"
            for it = 2:snt
                w = [-1.0 * srcwav[ipw][iss].d[sfield][it, is] for is = 1:ns]
                wavelets[nt-it+2][sfield] = Data.Array(w)
            end
        elseif (sflags[ipw] == 2)
            "ϕ[t] = ∑₁ᵗ⁻¹ s[t]"
            w = zeros(ns)
            if (sfield == :p)
                for it = 1:snt-1
                    for is = 1:ns
                        w[is] += (srcwav[ipw][iss].d[sfield][it, is] .* δt)
                    end
                    wavelets[it+1][sfield] = Data.Array(w)
                end
            else
                for it = 2:snt-1
                    for is = 1:ns
                        w[is] += (
                            (
                                (srcwav[ipw][iss].d[sfield][it, is] .* δt) +
                                (srcwav[ipw][iss].d[sfield][it-1, is] .* δt)
                            ) * 0.5
                        )
                    end
                    wavelets[it+1][sfield] = Data.Array(w)
                end

            end
            if (nt > snt)
                copyto!(wavelets[snt+1:end][sfield], wavelets[snt][sfield])
            end
        elseif (sflags[ipw] == -2)
            "source sink for 2"
            # multiplication with -1 for subtraction #time reversal
            # as the source wavelet has to be subtracted before the propagation step, I shift here by one sample"
            w = zeros(ns)
            for it = 1:snt-1
                for is = 1:ns
                    w[is] += (srcwav[ipw][iss].d[sfield][it, is] .* δt)
                end
                wavelets[nt-it+1][sfield][is] = Data.Array(-1.0 .* w)
            end
            if (nt > snt)
                nt_diff = nt - snt
                copyto!(wavelets[1:nt_diff+1][sfield], wavelets[nt_diff+2][sfield])
            end
        end
        w = vcat([Array(wavelets[it][sfield]) for it = 1:nt]...)
        freqmax =
            min(GeoPhyInv.Utils.findfreq(w, srcwav[ipw][iss].grid, attrib = :max), freqmax)
        freqmin =
            max(GeoPhyInv.Utils.findfreq(w, srcwav[ipw][iss].grid, attrib = :min), freqmin)
        push!(freqpeaks, GeoPhyInv.Utils.findfreq(w, srcwav[ipw][iss].grid, attrib = :peak))
    end
    return freqmin, freqmax, Statistics.mean(freqpeaks)
end

# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function add_source!(it::Int64, issp::Int64, iss::Int64, pac, pap)
    # adding source to respective sfield at [it] 
    for ipw in pac.activepw
        if (pac.sflags[ipw] ≠ 0) # add only if sflags is non-zero
            for sfield in pac.sfields[ipw]
                pw = view(pap[ipw].w1[:t][sfield], :)
                w = pap[ipw].ss[issp].wavelets[it][sfield]
                ssprayw = pap[ipw].ss[issp].ssprayw[sfield]
                mul!(pw, ssprayw, w, 1.0, 1.0)
            end
        end
    end
end


