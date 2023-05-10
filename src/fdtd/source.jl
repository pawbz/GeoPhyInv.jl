
# nothing # just put zeros, no sources added, when flag is put to zero
# -ve values correspond to source sink, used while time reversal
# here, size of w is (nt, ns)
get_source(w, ::Any, ::Val{0}) = zero(w)


get_source(w, ::p, ::Val{1}) = w
get_source(w, ::p, ::Val{-1}) = reverse!(w, dims = 1)
function get_source(w, ::Union{vz,vx,vy}, ::Val{1})
    ww = zero(w)
    for is = 1:size(ww, 2)
        for it = 1:size(ww, 1)
            # ww[it, is] = (w[it-1, is] + w[it, is]) * 0.5
            ww[it, is] =  w[it, is]
        end
    end
    return ww
end
function get_source(w, ::Union{vz,vx,vy}, ::Val{-1})
    ww = zero(w)
    for is = 1:size(ww, 2)
        for it = 1:size(ww, 1)
            # ww[it, is] = (w[it-1, is] + w[it, is]) * 0.5
            ww[it, is] = w[it, is]
        end
    end
    # as the source wavelet has to be subtracted before the propagation step, I shift here by one sample
    ww = circshift(ww, (-1, 0))
    reverse!(ww, dims = 1)
    ww[1, :] .= zero(Data.Number)
    # multiplication with -1 for subtraction #time reversal
    rmul!(ww, -one(Data.Number))
    return ww
end

# cumsum is used for integration, but `dt` factor is ignored
get_source(w, f::Any, ::Val{2}) = get_source(cumsum(w, dims = 1), f, Val{1}())
get_source(w, f::Any, ::Val{-2}) = get_source(cumsum(w, dims = 1), f, Val{-1}())



"""
Use input srcwav to fill wavelets and returns frequency bounds.
"""
function fill_wavelets!(ipw, iss, wavelets, srcwav, src_types)

    nt = size(wavelets, 1)
    sfields = names(srcwav[ipw][iss].d)[1]
    ns = srcwav[ipw][iss][:ns]
    δt = step(srcwav[ipw][iss].grid)

    freqmin = 0.0
    freqpeaks = []
    freqmax = Inf

    # change fieldnames that might be different
    for w in wavelets # time slice
        setnames!(w, sfields, 1)
    end

    for sfield in sfields
        w = [srcwav[ipw][iss].d[sfield][it, is] for it = 1:nt, is = 1:ns]
        w = get_source(w, eval(sfield)(), Val{src_types[ipw]}())
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
@inbounds @fastmath function add_source!(it::Int64, issp::Int64, iss::Int64, pac, pap, activepw, src_flags, sfields=Fields(pac.attrib_mod))
    # adding source to respective sfield at [it] 
    for ipw in activepw
        sfields = intersect(names(pac.srcwav[ipw][iss].d)[1], sfields)
        if (src_flags[ipw]) # add only if src_flags is non-zero
            for sfield in sfields
                pw = view(pap[ipw].w1[:t][sfield], :)
                w = pap[ipw].ss[issp].wavelets[it][sfield]
                ssprayw = pap[ipw].ss[issp].ssprayw[sfield]
                mul!(pw, ssprayw, w, 1.0, 1.0)
            end
        end
    end
end




"""
```julia
update!(pa,srcwav_new,src_types)
```
Update `pa` with a new bundle of source wavelets `srcwav_new`, without additional memory allocation.
Optionally, `src_types` can be changed. 

* `src_types=2` : source related flags 
  * `=0` inactive sources
  * `=1` sources with injection rate
  * `=2` volume injection sources
"""
function update!(pa::PFdtd, srcwav::SrcWav, src_types = 2; verbose = false)
    update!(pa, [srcwav], src_types; verbose = verbose)
end
function update!(pa::PFdtd, srcwav::Vector{SrcWav}, src_types = fill(1, length(srcwav)); verbose = false)
    # make a copy of srcwav in pa.c
    @assert (length(srcwav) == length(pa.c.srcwav))
    for i = 1:length(srcwav)
        copyto!(pa.c.srcwav[i], srcwav[i])
    end
    for ipw = 1:pa.c.ic[:npw]
        freqmin = 0.0
        freqmax = Inf
        freqpeaks = []
        # fill_wavelets for each supersource
        @sync begin
            for (ip, p) in enumerate(procs(pa.p))
                @async remotecall_wait(p) do
                    pap = localpart(pa.p)[ipw]
                    for is = 1:length(pap.ss)
                        iss = pap.ss[is].iss
                        wavelets = pap.ss[is].wavelets
                        broadcast(x -> fill!.(x, 0.0), wavelets)
                        fmin, fmax, fpeak =
                            fill_wavelets!(ipw, iss, wavelets, pa.c.srcwav, src_types)
                        freqmin = max(freqmin, fmin)
                        freqmax = min(freqmax, fmax)
                        push!(freqpeaks, fpeak)
                    end
                end
            end
        end

        # after srcwav is updated, sfields might have changed, so need to update spray and interpolation matrices as well
        update!(pa, pa.c.ageom, Srcs())

        freqpeak = Statistics.mean(freqpeaks)
        # store frequency (in Hz) bounds for first propagating wavefield
        if (ipw == 1)
            for (f, nf) in
                zip([freqmin, freqmax, freqpeak], [:freqmin, :freqmax, :freqpeak])
                if (nf ∈ names(pa.c.fc)[1])
                    pa.c.fc[nf] = Data.Number(f)
                else
                    pa.c.fc = vcat(pa.c.fc, NamedArray([Data.Number(f)], [nf]))
                end
            end
        end
        if (verbose)
            freqmin = @sprintf("%0.2e", freqmin)
            freqmax = @sprintf("%0.2e", freqmax)
            freqpeak = @sprintf("%0.2e", freqpeak)
            @info "frequency bounds for propagating wavefield $ipw are: [$freqmin, $freqmax], with peak at: $freqpeak"
        end
    end
end
