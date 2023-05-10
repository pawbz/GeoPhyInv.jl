for dimnames in [zip([:1, :2, :3], [:z, :y, :x]), zip([:1, :2], [:z, :x])]

    is = broadcast(x -> Symbol(string("i", x)), getindex.(collect(dimnames), 2))
    ist = Meta.parse(string("(", [string(s, ",") for s in is]..., ")"))
    N = Meta.parse(string(length(is)))

    for (idim, dim) in dimnames
        i = Symbol("i", string(dim))
        isboff = replace(is, i => :($i + boff))
        isdoff = replace(is, i => :(doff + $i))
        dimmin = Symbol(string(dim), "min")
        dimmax = Symbol(string(dim), "max")

        fnamehalf = Symbol("boundary_half", string(dim), "!")
        # d is the field
        # b is the boundary array
        @eval @parallel_indices($ist, function $fnamehalf(d::Data.Array{$N}, b, doff, boff)
            d[$(isdoff...)] = b[$(isboff...)]
            return
        end)

        fname = Symbol("boundary_force", string(dim), "!")
        @eval function $fname(d::Data.Array{$N}, b, pml_faces)
            np = ($(Meta.quot(dimmin)) ∈ pml_faces) ? _fd_npml : 0
            sb = collect(size(b))
            setindex!(sb, _fd_nbound, $idim)
            # first nbound points
            @parallel map(x -> (:)(1, x), Tuple(sb)) $fnamehalf(d, b, np, 0)
            # last nbound points independent of d
            @parallel map(x -> (:)(1, x), Tuple(sb)) $fnamehalf(
                d,
                b,
                getindex(size(d), $idim) - np - _fd_nbound,
                _fd_nbound,
            )
        end

        fname = Symbol("boundary_save", string(dim), "!")
        @eval function $fname(b::Data.Array{$N}, d, pml_faces)
            np = ($(Meta.quot(dimmin)) ∈ pml_faces) ? _fd_npml : 0
            sb = collect(size(b))
            setindex!(sb, _fd_nbound, $idim)
            # first _fd_np points
            @parallel map(x -> (:)(1, x), Tuple(sb)) $fnamehalf(b, d, 0, np)
            # last _fd_np points independent of d
            @parallel map(x -> (:)(1, x), Tuple(sb)) $fnamehalf(
                b,
                d,
                _fd_nbound,
                getindex(size(d), $idim) - np - _fd_nbound,
            )
        end

    end
end


# For each field, these functions output boundary storage arrays and a single snapshot.
for dimnames in [zip([:1, :2, :3], [:z, :y, :x]), zip([:1, :2], [:z, :x])]
    qs = broadcast(x -> Symbol(string(":", x)), getindex.(collect(dimnames), 2))
    dimsymbols = Meta.parse(string("[", [string(s, ",") for s in qs]..., ":snap]"))
    sizes1 = broadcast(x -> Symbol(string("n1", x)), getindex.(collect(dimnames), 2))
    sizes = broadcast(x -> Symbol(string("n", x)), getindex.(collect(dimnames), 2))
    sizes1_nbound = []
    ntvec = []
    for (idim, dim) in dimnames
        i = Symbol("n1", string(dim))
        ss = replace(sizes1, i => :(2 * _fd_nbound))
        push!(ntvec, [:nt]) # add time dimension
        push!(sizes1_nbound, ss)
    end
    push!(ntvec, [:1])
    push!(sizes1_nbound, sizes1) # for a single snapshot (no time dimension here)
    lsizes1_nbound = Meta.parse(string(length(sizes1_nbound)))

    for f in Fields(ndims=length(collect(dimnames)))
        @eval function get_boundary_store(::$f, attrib_mod, $(sizes...), nt)
            A = Array{Vector{Data.Array{_fd_ndims}}}(undef, $lsizes1_nbound)
            # get new sizes
            $(
                (
                    quote
                        $n1 = length(
                            get_mgrid(
                                $f(),
                                attrib_mod,
                                $((
                                    quote
                                        range(0, length=$n, stop=1)
                                    end for n in sizes
                                )...),
                            )[$i],
                        )
                    end for (i, n1) in enumerate(sizes1)
                )...
            )
            $(
                (
                    quote
                        A[$j] =
                            [Data.Array(zeros(($(sizes[1]...)))) for it = 1:$(sizes[2]...)]
                    end for (j, sizes) in enumerate(zip(sizes1_nbound, ntvec))
                )...
            )
            return NamedArray(A, $dimsymbols)
        end
    end
end



# for options not specified below, don't do anything
function boundary_force!(args...)
end
function boundary_force!(
    ::Val{:adjoint},
    it::Int64,
    issp::Int64,
    pac::T,
    pap,
) where {T<:P_common{<:FdtdAcoustic,2}}
    boundary = pap.ss[issp].boundary
    w1t = pap.w1[:t]
    boundary_forcex!(w1t[:p], boundary[:p][:x][it], pac.pml_faces)
    boundary_forcez!(w1t[:p], boundary[:p][:z][it], pac.pml_faces)
end
function boundary_force!(
    ::Val{:adjoint},
    it::Int64,
    issp::Int64,
    pac::T,
    pap,
) where {T<:P_common{<:FdtdElastic,2}}
    boundary = pap.ss[issp].boundary
    w1t = pap.w1[:t]
    for f in [:tauxx, :tauxz, :tauzz]
        boundary_forcex!(w1t[f], boundary[f][:x][it], pac.pml_faces)
        boundary_forcez!(w1t[f], boundary[f][:z][it], pac.pml_faces)
    end
end

function boundary_force!(
    ::Val{:adjoint},
    it::Int64,
    issp::Int64,
    pac::T,
    pap,
) where {T<:P_common{<:FdtdAcoustic,3}}
    boundary = pap.ss[issp].boundary
    w1t = pap.w1[:t]
    boundary_forcex!(w1t[:p], boundary[:p][:x][it], pac.pml_faces)
    boundary_forcey!(w1t[:p], boundary[:p][:y][it], pac.pml_faces)
    boundary_forcez!(w1t[:p], boundary[:p][:z][it], pac.pml_faces)
end
function boundary_force!(
    ::Val{:adjoint},
    it::Int64,
    issp::Int64,
    pac::T,
    pap,
) where {T<:P_common{<:FdtdElastic,3}}
    boundary = pap.ss[issp].boundary
    w1t = pap.w1[:t]
    for f in [:tauxx, :tauxz, :tauzz]
        boundary_forcex!(w1t[f], boundary[f][:x][it], pac.pml_faces)
        boundary_forcey!(w1t[f], boundary[f][:y][it], pac.pml_faces)
        boundary_forcez!(w1t[f], boundary[f][:z][it], pac.pml_faces)
    end
end


function boundary_force_snap_tau!(args...)
end
function boundary_force_snap_v!(args...)
end
function boundary_force_snap_tau!(::Val{:adjoint}, issp::Int64, pac::T, pap) where {T<:P_common{<:FdtdElastic}}
    w1t = pap.w1[:t]
    boundary = pap.ss[issp].boundary
    for f in [:tauxx, :tauxz, :tauzz]
        copyto!(w1t[f], boundary[f][:snap][1])
    end
end
function boundary_force_snap_tau!(::Val{:adjoint},
    issp::Int64,
    pac::T,
    pap,
) where {T<:P_common{<:FdtdAcoustic}}
    w1t = pap.w1[:t]
    boundary = pap.ss[issp].boundary
    copyto!(w1t[:p], boundary[:p][:snap][1])
end
function boundary_force_snap_v!(::Val{:adjoint},
    issp::Int64,
    pac::T,
    pap,
) where {T<:Union{P_common{<:FdtdAcoustic,2},P_common{<:FdtdElastic,2}}}
    w1t = pap.w1[:t]
    boundary = pap.ss[issp].boundary
    copyto!(w1t[:vx], boundary[:vx][:snap][1])
    copyto!(w1t[:vz], boundary[:vz][:snap][1])
end
function boundary_force_snap_v!(::Val{:adjoint},
    issp::Int64,
    pac::T,
    pap,
) where {T<:Union{P_common{<:FdtdAcoustic,3},P_common{<:FdtdElastic,3}}}
    w1t = pap.w1[:t]
    boundary = pap.ss[issp].boundary
    for f in [:vx, :vy, :vz]
        copyto!(w1t[f], boundary[f][:snap][1])
    end
end

# for options not specified below, don't do anything
function boundary_save!(args...)
end
function boundary_save!(
    ::Val{:forward},
    it::Int64,
    issp::Int64,
    pac::T,
    pap,
) where {T<:P_common{<:FdtdAcoustic,2}}
    boundary = pap.ss[issp].boundary
    w1t = pap.w1[:t]
    boundary_savex!(boundary[:p][:x][it], w1t[:p], pac.pml_faces)
    boundary_savez!(boundary[:p][:z][it], w1t[:p], pac.pml_faces)
    rmul!(boundary[:p][:x][it], -one(Data.Number))
    rmul!(boundary[:p][:z][it], -one(Data.Number))
end

function boundary_save!(
    ::Val{:forward},
    it::Int64,
    issp::Int64,
    pac::T,
    pap,
) where {T<:P_common{<:FdtdElastic,2}}
    boundary = pap.ss[issp].boundary
    w1t = pap.w1[:t]
    for f in [:tauxx, :tauxz, :tauzz]
        boundary_savex!(boundary[f][:x][it], w1t[f], pac.pml_faces)
        rmul!(boundary[f][:x][it], -one(Data.Number))
        boundary_savez!(boundary[f][:z][it], w1t[f], pac.pml_faces)
        rmul!(boundary[f][:z][it], -one(Data.Number))
    end
end

function boundary_save!(
    ::Val{:forward},
    it::Int64,
    issp::Int64,
    pac::T,
    pap,
) where {T<:P_common{<:FdtdAcoustic,3}}
    boundary = pap.ss[issp].boundary
    w1t = pap.w1[:t]
    boundary_savex!(boundary[:p][:x][it], w1t[:p], pac.pml_faces)
    boundary_savey!(boundary[:p][:y][it], w1t[:p], pac.pml_faces)
    boundary_savez!(boundary[:p][:z][it], w1t[:p], pac.pml_faces)
    rmul!(boundary[:p][:x][it], -one(Data.Number))
    rmul!(boundary[:p][:y][it], -one(Data.Number))
    rmul!(boundary[:p][:z][it], -one(Data.Number))
end


function boundary_save_snap_v!(args...)
end
function boundary_save_snap_tau!(args...)
end
function boundary_save_snap_tau!(::Val{:forward}, issp::Int64, pac::T, pap) where {T<:P_common{<:FdtdAcoustic}}
    w1t = pap.w1[:t]
    boundary = pap.ss[issp].boundary
    copyto!(boundary[:p][:snap][1], w1t[:p])
    rmul!(boundary[:p][:snap][1], -one(Data.Number))
end
function boundary_save_snap_tau!(::Val{:forward}, issp::Int64, pac::T, pap) where {T<:P_common{<:FdtdElastic}}
    w1t = pap.w1[:t]
    boundary = pap.ss[issp].boundary
    for f in [:tauxx, :tauxz, :tauzz]
        copyto!(boundary[f][:snap][1], w1t[f])
        rmul!(boundary[f][:snap][1], -one(Data.Number))
    end
end
function boundary_save_snap_v!(::Val{:forward},
    issp::Int64,
    pac::T,
    pap,
) where {T<:Union{P_common{<:FdtdAcoustic,2},P_common{<:FdtdElastic,2}}}
    w1t = pap.w1[:t]
    boundary = pap.ss[issp].boundary
    for f in [:vx, :vz]
        copyto!(boundary[f][:snap][1], w1t[f])
    end
end
function boundary_save_snap_v!(::Val{:forward},
    issp::Int64,
    pac::T,
    pap,
) where {T<:Union{P_common{<:FdtdAcoustic,3},P_common{<:FdtdElastic,3}}}
    w1t = pap.w1[:t]
    boundary = pap.ss[issp].boundary
    for f in [:vx, :vy, :vz]
        copyto!(boundary[f][:snap][1], w1t[f])
    end
end

