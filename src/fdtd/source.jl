
# nothing # just put zeros, no sources added, when flag is put to zero
# -ve values correspond to source sink, used while time reversal
# here, size of w is (nt, ns)
get_source(w, ::Any, ::Val{0}) = zero(w)


get_source(w, ::p, ::Val{1}) = w
get_source(w, ::p, ::Val{-1}) = reverse!(w, dims = 1)
get_source(w, ::p, ::Val{-1}) = reverse!(w, dims = 1)
function get_source(w, ::Union{vz,vx,vy}, ::Val{1})
    ww = zero(w)
    for is = 1:size(ww, 2)
        for it = 2:size(ww, 1)
            ww[it, is] = (w[it-1, is] + w[it, is]) * 0.5
        end
    end
    return ww
end
function get_source(w, ::Union{vz,vx,vy}, ::Val{-1})
    ww = zero(w)
    for is = 1:size(ww, 2)
        for it = 2:size(ww, 1)
            # as the source wavelet has to be subtracted before the propagation step, I shift here by one sample
            ww[it, is] = (w[it-1, is] + w[it, is]) * 0.5
        end
    end
    ww = circshift(ww, (1, 0))
    ww[1, :] .= zero(Data.Number)
    reverse!(ww, dims = 1)
    # multiplication with -1 for subtraction #time reversal
    rmul!(ww, -one(Data.Number))
    return ww
end

get_source(w, f::Any, ::Val{2}) = get_source(cumsum(w, dims = 1), f, Val{1}())
get_source(w, f::Any, ::Val{-2}) = get_source(cumsum(w, dims = 1), f, Val{-1}())



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
        w = [srcwav[ipw][iss].d[sfield][it, is] for it = 1:nt, is = 1:ns]
        w = get_source(w, eval(sfield)(), Val{sflags[ipw]}())
        snt = length(srcwav[ipw][1].grid)
        for it = 1:snt
            ww = view(w, it, :) # view for all sources
            copyto!(wavelets[it][sfield], Data.Array(ww))
        end
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


