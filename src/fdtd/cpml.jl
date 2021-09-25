"""
Generate coefficients related to PML boundaries, e.g., damping profiles.
This function outputs a, b, and k variables in eq. 25-26 from Komatitsch 2007 (Geophysics).

	TODO
* At the moment, K=1, and is dummy, unlike in the case of EM when it was useful. Removing K and K_inv might save us memory when using GPUs.
"""
function update_pml!(
    pml::NamedVector{T},
    mgrid, # 1D grid
    flags::Vector{Bool},
    δt::Float64,
    np::Int64,
    velavg::Float64,
    freqpeak::Float64, # dominant frequency in Hz
) where {T<:Data.Array}

    δx = step(mgrid)
    nx = length(mgrid)

    NPOWER = 2.e0
    K_MAX_PML = 1.e0 # from Gedney page 8.11
    ALPHA_MAX_PML = pi * freqpeak # from Festa and Vilotte

    # thickness of the PML layer in meters
    thickness_PML = np * δx

    "reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf"
    Rcoef = 0.001e0

    #! check that NPOWER is okay
    #if(NPOWER < 1) stop "NPOWER must be greater than 1"

    # compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
    d0 = -(NPOWER + 1) * (velavg) * log(Rcoef) / (2.e0 * thickness_PML)

    k = zeros(nx)
    k_half = zero(k)
    fill!(k, 1)
    fill!(k_half, 1)
    d = zeros(nx)
    d_half = zero(d)
    alpha = zero(d)
    alpha_half = zero(d)
    a=zero(d)
    b=zero(d)
    a_half=zero(d)
    b_half=zero(d)


    "damping in the X direction"
    "origin of the PML layer (position of right edge minus thickness, in meters)"
    xoriginleft = thickness_PML
    xoriginright = (nx - 1) * δx - thickness_PML

    for ix = 1:nx
        # abscissa of current grid point along the damping profile
        xval = δx * real(ix - 1)

        #---------- left edge
        if (flags[1])

            # define damping profile at the grid points
            abscissa_in_PML = xoriginleft - xval
            if (abscissa_in_PML >= 0.0)
                abscissa_normalized = abscissa_in_PML / thickness_PML
                d[ix] = d0 * abscissa_normalized .^ NPOWER
                # this taken from Gedney page 8.2
                k[ix] = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized .^ NPOWER
                alpha[ix] =
                    ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.01e0 * ALPHA_MAX_PML
            end

            # define damping profile at half the grid points
            abscissa_in_PML = xoriginleft - (xval + δx / 2.e0)
            if (abscissa_in_PML >= 0.0)
                abscissa_normalized = abscissa_in_PML / thickness_PML
                d_half[ix] = d0 * abscissa_normalized .^ NPOWER
                # this taken from Gedney page 8.2
                k_half[ix] = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized .^ NPOWER
                alpha_half[ix] =
                    ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.01e0 * ALPHA_MAX_PML
            end
        end

        #---------- right edge
        if (flags[2])

            # define damping profile at the grid points
            abscissa_in_PML = xval - xoriginright
            if (abscissa_in_PML >= 0.0)
                abscissa_normalized = abscissa_in_PML / thickness_PML
                d[ix] = d0 * abscissa_normalized^NPOWER
                # this taken from Gedney page 8.2
                k[ix] = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized^NPOWER
                alpha[ix] =
                    ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.01e0 * ALPHA_MAX_PML
            end

            # define damping profile at half the grid points
            abscissa_in_PML = xval + δx / 2.e0 - xoriginright
            if (abscissa_in_PML >= 0.0)
                abscissa_normalized = abscissa_in_PML / thickness_PML
                d_half[ix] = d0 * abscissa_normalized^NPOWER
                # this taken from Gedney page 8.2
                k_half[ix] = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized^NPOWER
                alpha_half[ix] =
                    ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.01e0 * ALPHA_MAX_PML
            end

        end
        #
        # just in case, for -5 at the end
        (alpha[ix] < 0.0) ? alpha[ix] = 0.0 : nothing
        (alpha_half[ix] < 0.0) ? alpha_half[ix] = 0.0 : nothing

        # see equation 25 Komatitsch, 2007, get b and a (need k, alpha, and d before this )
        b[ix] = exp(-(d[ix] / k[ix] + alpha[ix]) * δt)
        b_half[ix] = exp(-(d_half[ix] / k_half[ix] + alpha_half[ix]) * δt)

        # this to avoid division by zero outside the PML
        (abs(d[ix]) > 1.e-6) ?
        a[ix] = d[ix] * (b[ix] - 1.e0) / (k[ix] * (d[ix] + k[ix] * alpha[ix])) : nothing
        (abs(d_half[ix]) > 1.e-6) ?
        a_half[ix] =
            d_half[ix] * (b_half[ix] - 1.e0) /
            (k_half[ix] * (d_half[ix] + k_half[ix] * alpha_half[ix])) : nothing

    end
    copyto!(pml[:kI], inv.(k))
    copyto!(pml[:k_halfI], inv.(k_half))
    copyto!(pml[:a],a)
    copyto!(pml[:b],b)
    copyto!(pml[:a_half],a_half)
    copyto!(pml[:b_half],b_half)

end


"""
Generate a NamedArray with PML coefficients for all the dimensions that are then stored in the FDTD structs.
"""
function get_pml(mgrid)
    dnames = dim_names(length(mgrid))
    pnames = [:a, :b, :kI, :a_half, :b_half, :k_halfI]
    return NamedArray(
        [
            NamedArray(
                Data.Array.([
                    zeros(ni),
                    zeros(ni),
                    ones(ni),
                    zeros(ni),
                    zeros(ni),
                    ones(ni),
                ]),
                pnames,
            ) for ni in length.(mgrid)
        ],
        dnames,
    )
end

"""
Just a loop over dims of pml
"""
function update_pml!(
    pml::NamedVector{T},
    mgrid,
    pml_edges::Vector{Symbol},
    args...,
) where {T<:NamedVector}
    for (i, dim) in enumerate(dim_names(length(mgrid)))
        update_pml!(
            pml[dim],
            mgrid[i],
            [
                any(pml_edges .== Symbol(string(dim), "min")),
                any(pml_edges .== Symbol(string(dim), "max")),
            ],
            args...,
        )
    end
end

function update_pml!(pac)
    update_pml!(
        pac.pml,
        pac.exmedium.mgrid,
        pac.pml_edges,
        pac.fc[:dt],
        npml-3,
        Statistics.mean(pac.exmedium.bounds[:vp]),
        pac.fc[:freqpeak],
    )
end


for dimnames in [zip([:1, :2, :3], dim_names(3)), zip([:1, :2], dim_names(2))]
    is = broadcast(x -> Symbol(string("i", x)), getindex.(collect(dimnames), 2))
    ist = Meta.parse(string("(", [string(s, ",") for s in is]..., ")"))
    N = Meta.parse(string(length(is)))

    for (idim, dim) in dimnames
        i = Symbol("i", string(dim))
        ismoff = replace(is, i => :($i + moff))
        isdoff = replace(is, i => :(doff + $i))
        for (fname, fnamenp, idoff) in zip(
            [Symbol("memory", string(dim), "!"), Symbol("memory1", string(dim), "!")],
            [Symbol("memorynp", string(dim), "!"), Symbol("memorynp1", string(dim), "!")],
            [:($i + doff), :($i + doff + 1)],
        )
            @eval @parallel_indices(
                $ist,
                function $fnamenp(memory::Data.Array{$N}, d, a, b, kI, moff, doff)
                    memory[$(ismoff...)] =
                        b[$idoff] * memory[$(ismoff...)] + a[$idoff] * d[$(isdoff...)]
                    d[$(isdoff...)] = d[$(isdoff...)] * kI[$idoff] + memory[$(ismoff...)]
                    return
                end
            )
            @eval function $fname(memory::Data.Array{$N}, d, a, b, kI)
                sm = collect(size(memory))
                setindex!(sm, npml, $idim)
                # first npml points
                @parallel map(x -> (:)(1, x), Tuple(sm)) $fnamenp(memory, d, a, b, kI, 0, 0)
                # last npml points independent of d
                @parallel map(x -> (:)(1, x), Tuple(sm)) $fnamenp(
                    memory,
                    d,
                    a,
                    b,
                    kI,
                    npml,
                    getindex(size(d), $idim) - npml,
                )
            end
        end
    end
end

