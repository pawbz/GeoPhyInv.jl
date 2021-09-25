"""
Use input srcwav to fill wavelets and returns frequency bounds.
"""
function fill_wavelets!(ipw, iss, wavelets, srcwav, sflags)

    nt = size(wavelets, 1)
    sfields = names(srcwav[ipw][iss].d)[1]
    ns = srcwav[ipw][iss][:ns]
    δt = step(srcwav[ipw][iss].grid)

    freqmin = 0.0
    freqpeaks= []
    freqmax = Inf
    for sfield in sfields
        for is = 1:ns
            snt = length(srcwav[ipw][1].grid)
            if (sflags[ipw] == 0)
                # nothing # just put zeros, no sources added
                for it = 1:nt
                    wavelets[it][sfield][is] = 0.0
                end
            elseif (sflags[ipw] == 1)
                "ϕ[t] = s[t]"
                for it = 1:snt
                    source_term = srcwav[ipw][iss].d[sfield][it, is]
                    wavelets[it][sfield][is] = source_term
                end
            elseif (sflags[ipw] == -1)
                "source sink for 1"
                "ϕ[t] = -s[t]"
                for it = 2:snt
                    source_term = srcwav[ipw][iss].d[sfield][it, is]
                    wavelets[nt-it+2][sfield][is] = -1.0 * source_term
                end
            elseif (sflags[ipw] == 2)
                "ϕ[t] = ∑₁ᵗ⁻¹ s[t]"
                source_term_stack = 0.0
                if (sfield == :p)
                    for it = 1:snt-1
                        source_term_stack += (srcwav[ipw][iss].d[sfield][it, is] .* δt)
                        wavelets[it+1][sfield][is] = source_term_stack
                    end
                else
                    for it = 2:snt-1
                        source_term_stack += (
                            (
                                (srcwav[ipw][iss].d[sfield][it, is] .* δt) +
                                (srcwav[ipw][iss].d[sfield][it-1, is] .* δt)
                            ) * 0.5
                        )
                        wavelets[it+1][sfield][is] = source_term_stack
                    end

                end
                if (nt > snt)
                    wavelets[snt+1:end][sfield][is] = wavelets[snt][sfield][is]
                end
            elseif (sflags[ipw] == -2)
                "source sink for 2"
                # multiplication with -1 for subtraction #time reversal
                # as the source wavelet has to be subtracted before the propagation step, I shift here by one sample"
                source_term_stack = 0.0
                for it = 1:snt-1
                    source_term_stack += (srcwav[ipw][iss].d[sfield][it, is] .* δt)
                    wavelets[nt-it+1][sfield][is] = -1.0 * source_term_stack
                end
                if (nt > snt)
                    nt_diff = nt - snt
                    wavelets[1:nt_diff+1][sfield][is] = wavelets[nt_diff+2][sfield][is]
                end
            end
            w = [wavelets[it][sfield][is] for it = 1:nt]
            freqmax = min(
                GeoPhyInv.Utils.findfreq(w, srcwav[ipw][iss].grid, attrib = :max),
                freqmax,
            )
            freqmin = max(
                GeoPhyInv.Utils.findfreq(w, srcwav[ipw][iss].grid, attrib = :min),
                freqmin,
            )
            push!(freqpeaks,GeoPhyInv.Utils.findfreq(w, srcwav[ipw][iss].grid, attrib = :peak))
        end
    end
    return freqmin, freqmax, Statistics.mean(freqpeaks)
end

struct Source_B1 end
struct Source_B0 end

# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function add_source!(it::Int64, issp::Int64, iss::Int64, pac, pap)
    # aliases
    ageom = pac.ageom
    """
    adding source to pressure field at [it] 
    """
    for ipw in pac.activepw
        p = pap[ipw].w1[:t]
        wavelets = pap[ipw].ss[issp].wavelets
        sindices = pap[ipw].ss[issp].sindices
        ssprayw = pap[ipw].ss[issp].ssprayw
        sfields = pac.sfields[ipw]
        if (pac.sflags[ipw] ≠ 0) # only if sflags is non-zero
            for sfield in sfields
                @simd for is = 1:ageom[ipw][iss].ns
                    """
                    use wavelets at [it], i.e., sum of source terms
                    until [it-1]
                    division of source term with δx and δz (see Jan's fdelmodc manual)
                    """
                    source_term =
                        wavelets[it][sfield][is] *
                        pac.fc[:dt] *
                        prod(inv.(step.(pac.medium.mgrid)))

                    for (i, si) in enumerate(sindices[is])
                        p[sfield][si] +=
                            source(source_term, ssprayw[is][i], pac, si, eval(sfield)())
                    end
                end
            end
        end
    end
end


# multiplication with modK
# on pressure grid
source(source_term, spray, pac::P_common{FdtdAcou}, si, ::p) =
    source_term * spray * pac.mod[:K][si]

# on vx grid
source(source_term, spray, pac::P_common{FdtdAcou}, si, ::vx) = source_term * spray #* pac.mod[:rhovxI][si]
source(source_term, spray, pac::P_common{FdtdElastic}, si, ::vx) = source_term * spray #* pac.mod[:rhovxI][si]

# on vz grid
source(source_term, spray, pac::P_common{FdtdAcou}, si, ::vz) = source_term * spray #* pac.mod[:rhovzI][si]
source(source_term, spray, pac::P_common{FdtdElastic}, si, ::vz) = source_term * spray #* pac.mod[:rhovzI][si]



# testing source on tauxx (need to check it)
source(source_term, spray, pac::P_common{FdtdElastic}, si, ::tauxx) =
    source_term * spray * 10^14
source(source_term, spray, pac::P_common{FdtdElastic}, si, ::tauyy) =
    source_term * spray * 10^14
