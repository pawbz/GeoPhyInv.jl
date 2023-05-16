"""
Generate coefficients related to PML boundaries, e.g., damping profiles.
This function outputs a, b, and k variables in eq. 25-26 from Komatitsch 2007 (Geophysics).

	TODO
* At the moment, K=1, and is dummy, unlike in the case of EM when it was useful. Removing K and K_inv might save us memory when using GPUs.
"""
function update_pml!(
    pml::NamedVector{T},
    exmgrid, # 1D grid extended
    mgrid, # 1D grid without extension
    flags::Vector{Bool},
    δt,
    velavg,
    freqpeak, # dominant frequency in Hz
) where {T<:Data.Array}

    nx = length(exmgrid)

    # origin, where the PMLs start (offset b/w the domain of interest and PML region depends on order
    xoriginleft = mgrid[1] - 2 * (step(mgrid))
    xoriginright = mgrid[end] + 2 * (step(mgrid))

    NPOWER = 2.e0
    K_MAX_PML = 1.e0 # from Gedney page 8.11
    ALPHA_MAX_PML = pi * freqpeak # from Festa and Vilotte

    # thickness of the PML layer in meters
    thickness_PML = (_fd_npml - 1) * step(mgrid)

    "reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf"
    Rcoef = 0.001e0

    #! check that NPOWER is okay
    #if(NPOWER < 1) stop "NPOWER must be greater than 1"

    # compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
    d0 = -(NPOWER + 1) * (velavg) * log(Rcoef) / (2.e0 * thickness_PML)

    k = zeros(nx)
    fill!(k, 1)
    d = zeros(nx)
    alpha = zero(d)
    a = zero(d)
    b = zero(d)

    "damping in the X direction"
    "origin of the PML layer (position of right edge minus thickness, in meters)"

    for ix = 1:nx
        #---------- left edge
        if (flags[1])
            # define damping profile at the grid points
            abscissa_in_PML = xoriginleft - exmgrid[ix]
            if (abscissa_in_PML >= 0.0)
                abscissa_normalized = abscissa_in_PML / thickness_PML
                d[ix] = d0 * abscissa_normalized .^ NPOWER
                # this taken from Gedney page 8.2
                k[ix] = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized .^ NPOWER
                alpha[ix] =
                    ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.01e0 * ALPHA_MAX_PML
            end
        end
        #---------- right edge
        if (flags[2])
            # define damping profile at the grid points
            abscissa_in_PML = exmgrid[ix] - xoriginright
            if (abscissa_in_PML >= 0.0)
                abscissa_normalized = abscissa_in_PML / thickness_PML
                d[ix] = d0 * abscissa_normalized^NPOWER
                # this taken from Gedney page 8.2
                k[ix] = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized^NPOWER
                alpha[ix] =
                    ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.01e0 * ALPHA_MAX_PML
            end
        end
        #
        # just in case, for -5 at the end
        (alpha[ix] < 0.0) ? alpha[ix] = 0.0 : nothing

        # see equation 25 Komatitsch, 2007, get b and a (need k, alpha, and d before this )
        b[ix] = exp(-(d[ix] / k[ix] + alpha[ix]) * δt)

        # this to avoid division by zero outside the PML
        (abs(d[ix]) > 1.e-6) ?
        a[ix] = d[ix] * (b[ix] - 1.e0) / (k[ix] * (d[ix] + k[ix] * alpha[ix])) : nothing
    end
    pkI = zeros(2 * _fd_npml)
    pa = zeros(2 * _fd_npml)
    pb = zeros(2 * _fd_npml)
    for i = 1:_fd_npml
        j = i + _fd_npml
        pkI[i] = inv(k[i])
        pkI[j] = inv(k[end-_fd_npml+i])
        pa[i] = a[i]
        pa[j] = a[end-_fd_npml+i]
        pb[i] = b[i]
        pb[j] = b[end-_fd_npml+i]
    end
    copyto!(pml[:kI], pkI)
    copyto!(pml[:a], pa)
    copyto!(pml[:b], pb)
end


"""
Generate a NamedArray with PML coefficients for all the dimensions that are then stored in the FDTD structs.
"""
function get_pml(attrib_mod, mgrid)
    dfields = Fields(attrib_mod, "d", ndims = length(mgrid)) # derivative fields that need PML memory
    pnames = [:a, :b, :kI]
    np = 2 * _fd_npml
    return NamedArray(
        [
            NamedArray(Data.Array.([zeros(np), zeros(np), ones(np)]), pnames) for
            df in dfields
        ],
        dfields,
    )
end

"""
Just a loop over dims of pml
"""
function update_pml!(
    pml::NamedVector{T},
    exmgrid,
    mgrid,
    pml_faces::Vector{Symbol},
    attrib_mod,
    args...,
) where {T<:NamedVector}
    for df in names(pml)[1]
        dim = string(last(string(df)))
        i = findfirst(x -> string(x) == dim, dim_names(length(exmgrid)))
        exmgridf = get_mgrid(eval(df)(), attrib_mod, exmgrid...)[i]
        mgridf = mgrid[i]
        flags =
            [Symbol(string(dim), "min") ∈ pml_faces, Symbol(string(dim), "max") ∈ pml_faces]
        update_pml!(pml[df], exmgridf, mgridf, flags, args...)
    end
end

function update_pml!(pac)
    update_pml!(
        pac.pml,
        pac.exmedium.mgrid,
        pac.medium.mgrid,
        pac.pml_faces,
        pac.attrib_mod,
        pac.fc[:dt],
        Statistics.mean(pac.exmedium[:vp].bounds),
        pac.fc[:freqpeak],
    )
end


for dimnames in [zip([:1, :2, :3], dim_names(3)), zip([:1, :2], dim_names(2))]
    is = broadcast(x -> Symbol(string("i", x)), getindex.(collect(dimnames), 2))
    ist = Meta.parse(string("(", [string(s, ",") for s in is]..., ")"))
    N = Meta.parse(string(length(is)))

    for (idim, dim) in dimnames
        i = Symbol("i", string(dim))
        dimmin = Symbol(string(dim), "min")
        dimmax = Symbol(string(dim), "max")

        ismoff = replace(is, i => :($i + moff))
        isdoff = replace(is, i => :(doff + $i))
        for (fname, fnamenp, imoff) in zip(
            [Symbol("memory", string(dim), "!"), Symbol("memory1", string(dim), "!")],
            [Symbol("memorynp", string(dim), "!"), Symbol("memorynp1", string(dim), "!")],
            [:($i + moff), :($i + moff + 1)],
        )
            @eval @parallel_indices(
                $ist,
                function $fnamenp(memory::Data.Array{$N}, d, a, b, kI, moff, doff)
                    memory[$(ismoff...)] =
                        b[$imoff] * memory[$(ismoff...)] + a[$imoff] * d[$(isdoff...)]
                    d[$(isdoff...)] = d[$(isdoff...)] * kI[$imoff] + memory[$(ismoff...)]
                    return
                end
            )
            @eval function $fname(memory::Data.Array{$N}, d, a, b, kI, pml_faces)

                sm = collect(size(memory))
                setindex!(sm, _fd_npml, $idim)
                # first _fd_npml points
                if ($(Meta.quot(dimmin)) ∈ pml_faces)
                    @parallel map(x -> (:)(1, x), Tuple(sm)) $fnamenp(
                        memory,
                        d,
                        a,
                        b,
                        kI,
                        0,
                        0,
                    )
                end
                # last _fd_npml points independent of d
                if ($(Meta.quot(dimmax)) ∈ pml_faces)
                    @parallel map(x -> (:)(1, x), Tuple(sm)) $fnamenp(
                        memory,
                        d,
                        a,
                        b,
                        kI,
                        _fd_npml,
                        getindex(size(d), $idim) - _fd_npml,
                    )
                end
            end
        end
    end
end

