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
            np = ($(Meta.quot(dimmin)) ∈ pml_faces) ? _fd.npml : 0
            sb = collect(size(b))
            sb = collect(size(b))
            setindex!(sb, _fd.nbound, $idim)
            # first _fd.npml points
            @parallel map(x -> (:)(1, x), Tuple(sb)) $fnamehalf(d, b, np, 0)
            # last _fd.npml points independent of d
            @parallel map(x -> (:)(1, x), Tuple(sb)) $fnamehalf(
                d,
                b,
                getindex(size(d), $idim) - np - _fd.nbound,
                _fd.nbound,
            )
        end

        fname = Symbol("boundary_save", string(dim), "!")
        @eval function $fname(b::Data.Array{$N}, d, pml_faces)
            np = ($(Meta.quot(dimmin)) ∈ pml_faces) ? _fd.npml : 0
            sb = collect(size(b))
            setindex!(sb, _fd.nbound, $idim)
            # first _fd.npml points
            @parallel map(x -> (:)(1, x), Tuple(sb)) $fnamehalf(b, d, 0, np)
            # last _fd.npml points independent of d
            @parallel map(x -> (:)(1, x), Tuple(sb)) $fnamehalf(
                b,
                d,
                _fd.nbound,
                getindex(size(d), $idim) - np - _fd.nbound,
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
        ss = replace(sizes1, i => :(2 * _fd.nbound))
        push!(ntvec, [:nt]) # add time dimension
        push!(sizes1_nbound, ss)
    end
    push!(ntvec, [:1])
    push!(sizes1_nbound, sizes1) # for a single snapshot (no time dimension here)
    lsizes1_nbound = Meta.parse(string(length(sizes1_nbound)))

    for f in Fields(ndims = length(collect(dimnames)))
        @eval function get_boundary_store(::$f, attrib_mod, $(sizes...), nt)
            A = Array{Vector{Data.Array{_fd.ndims}}}(undef, $lsizes1_nbound)
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
                                            range(0, length = $n, stop = 1)
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


function boundary_force_snap_tau!(
    issp::Int64,
    pac::T,
    pap,
) where {T<:P_common{FdtdElastic}}
    w1t = pap.w1[:t]
    boundary = pap.ss[issp].boundary
    for f in [:tauxx, :tauxz, :tauzz]
        copyto!(w1t[f], boundary[f][:snap][1])
    end
end
function boundary_force_snap_tau!(
    issp::Int64,
    pac::T,
    pap,
) where {T<:P_common{FdtdAcoustic}}
    w1t = pap.w1[:t]
    boundary = pap.ss[issp].boundary
    copyto!(w1t[:p], boundary[:p][:snap][1])
end
function boundary_force_snap_v!(
    issp::Int64,
    pac::T,
    pap,
) where {T<:Union{P_common{FdtdAcoustic,2},P_common{FdtdElastic,2}}}
    w1t = pap.w1[:t]
    boundary = pap.ss[issp].boundary
    copyto!(w1t[:vx], boundary[:vx][:snap][1])
    copyto!(w1t[:vz], boundary[:vz][:snap][1])
end
function boundary_force_snap_v!(
    issp::Int64,
    pac::T,
    pap,
) where {T<:Union{P_common{FdtdAcoustic,3},P_common{FdtdElastic,3}}}
    w1t = pap.w1[:t]
    boundary = pap.ss[issp].boundary
    for f in [:vx, :vy, :vz]
    copyto!(w1t[f], boundary[f][:snap][1])
    end
end

function boundary_force!(
    it::Int64,
    issp::Int64,
    pac::T,
    pap,
) where {T<:P_common{FdtdAcoustic,2}}
    boundary = pap.ss[issp].boundary
    w1t = pap.w1[:t]
    boundary_forcex!(w1t[:p], boundary[:p][:x][it], pac.pml_faces)
    boundary_forcez!(w1t[:p], boundary[:p][:z][it], pac.pml_faces)
end
function boundary_force!(
    it::Int64,
    issp::Int64,
    pac::T,
    pap,
) where {T<:P_common{FdtdElastic,2}}
    boundary = pap.ss[issp].boundary
    w1t = pap.w1[:t]
    for f in [:tauxx, :tauxz, :tauzz]
    boundary_forcex!(w1t[f], boundary[f][:x][it], pac.pml_faces)
    boundary_forcez!(w1t[f], boundary[f][:z][it], pac.pml_faces)
    end
end

function boundary_force!(
    it::Int64,
    issp::Int64,
    pac::T,
    pap,
) where {T<:P_common{FdtdAcoustic,3}}
    boundary = pap.ss[issp].boundary
    w1t = pap.w1[:t]
    boundary_forcex!(w1t[:p], boundary[:p][:x][it], pac.pml_faces)
    boundary_forcey!(w1t[:p], boundary[:p][:y][it], pac.pml_faces)
    boundary_forcez!(w1t[:p], boundary[:p][:z][it], pac.pml_faces)
end
function boundary_force!(
    it::Int64,
    issp::Int64,
    pac::T,
    pap,
) where {T<:P_common{FdtdElastic,3}}
    boundary = pap.ss[issp].boundary
    w1t = pap.w1[:t]
    for f in [:tauxx, :tauxz, :tauzz]
        boundary_forcex!(w1t[f], boundary[f][:x][it], pac.pml_faces)
        boundary_forcey!(w1t[f], boundary[f][:y][it], pac.pml_faces)
        boundary_forcez!(w1t[f], boundary[f][:z][it], pac.pml_faces)
    end
end


function boundary_save_snap_tau!(issp::Int64, pac::T, pap) where {T<:P_common{FdtdAcoustic}}
    w1t = pap.w1[:t]
    boundary = pap.ss[issp].boundary
    copyto!(boundary[:p][:snap][1], w1t[:p])
    rmul!(boundary[:p][:snap][1], -one(Data.Number))
end
function boundary_save_snap_tau!(issp::Int64, pac::T, pap) where {T<:P_common{FdtdElastic}}
    w1t = pap.w1[:t]
    boundary = pap.ss[issp].boundary
    for f in [:tauxx, :tauxz, :tauzz]
        copyto!(boundary[f][:snap][1], w1t[f])
        rmul!(boundary[f][:snap][1], -one(Data.Number))
    end
end
function boundary_save_snap_v!(
    issp::Int64,
    pac::T,
    pap,
) where {T<:Union{P_common{FdtdAcoustic,2},P_common{FdtdElastic,2}}}
    w1t = pap.w1[:t]
    boundary = pap.ss[issp].boundary
    for f in [:vx, :vz]
        copyto!(boundary[f][:snap][1], w1t[f])
    end
end
function boundary_save_snap_v!(
    issp::Int64,
    pac::T,
    pap,
) where {T<:Union{P_common{FdtdAcoustic,3},P_common{FdtdElastic,3}}}
    w1t = pap.w1[:t]
    boundary = pap.ss[issp].boundary
    for f in [:vx, :vy, :vz]
        copyto!(boundary[f][:snap][1], w1t[f])
    end
end

function boundary_save!(
    it::Int64,
    issp::Int64,
    pac::T,
    pap,
) where {T<:P_common{FdtdAcoustic,2}}
    boundary = pap.ss[issp].boundary
    w1t = pap.w1[:t]
    boundary_savex!(boundary[:p][:x][it], w1t[:p], pac.pml_faces)
    boundary_savez!(boundary[:p][:z][it], w1t[:p], pac.pml_faces)
    rmul!(boundary[:p][:x][it], -one(Data.Number))
    rmul!(boundary[:p][:z][it], -one(Data.Number))
end

function boundary_save!(
    it::Int64,
    issp::Int64,
    pac::T,
    pap,
) where {T<:P_common{FdtdElastic,2}}
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
    it::Int64,
    issp::Int64,
    pac::T,
    pap,
) where {T<:P_common{FdtdAcoustic,3}}
    boundary = pap.ss[issp].boundary
    w1t = pap.w1[:t]
    boundary_savex!(boundary[:p][:x][it], w1t[:p], pac.pml_faces)
    boundary_savey!(boundary[:p][:y][it], w1t[:p], pac.pml_faces)
    boundary_savez!(boundary[:p][:z][it], w1t[:p], pac.pml_faces)
    rmul!(boundary[:p][:x][it], -one(Data.Number))
    rmul!(boundary[:p][:y][it], -one(Data.Number))
    rmul!(boundary[:p][:z][it], -one(Data.Number))
end



# @fastmath @inbounds function boundarys_l!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
# 	pw=p[:p]
# 	for ix=1:3
# 		@simd for iz=1:ibz1-ibz0+1
# 			boundary[4][iz,ix,it] = pw[ibz0+iz-1,ibx0+ix-1] 
# 		end
# 	end
# end
# @fastmath @inbounds function boundarys_r!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
# 	pw=p[:p]
# 	for ix=1:3
# 		@simd for iz=1:ibz1-ibz0+1
# 			boundary[2][iz,ix,it] = pw[ibz0+iz-1,ibx1-ix+1]
# 		end
# 	end
# end
# @fastmath @inbounds function boundarys_t!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
# 	pw=p[:p]
# 	@simd for ix=1:ibx1-ibx0+1
# 		for iz=1:3
# 			boundary[1][iz,ix,it] = pw[ibz0+iz-1,ibx0+ix-1]
# 		end
# 	end
# end
# @fastmath @inbounds function boundarys_b!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
# 	pw=p[:p]
# 	@simd for ix=1:ibx1-ibx0+1
# 		for iz=1:3
# 			boundary[3][iz,ix,it] = pw[ibz1-iz+1,ibx0+ix-1]
# 		end
# 	end
# end


# @fastmath @inbounds function boundary_force!(it::Int64,issp::Int64,pac,pap)
# 	boundary=pap[1].ss[issp].boundary
# 	p=pap[1].w1[:t]
# 	ibx0=pac.bindices[:bx0]; ibz0=pac.bindices[:bz0]; ibx1=pac.bindices[:bx1]; ibz1=pac.bindices[:bz1]
# 	boundaryf_l!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
# 	boundaryf_r!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
# 	boundaryf_t!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
# 	boundaryf_b!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
# end
# @fastmath @inbounds function boundaryf_l!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
# 	pw=p[:p]
# 	for ix=1:3
# 		@simd for iz=1:ibz1-ibz0+1
# 			pw[ibz0+iz-1,ibx0+ix-1] = boundary[4][iz,ix,it]
# 		end
# 	end
# end
# @fastmath @inbounds function boundaryf_r!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
# 	pw=p[:p]
# 	for ix=1:3
# 		@simd for iz=1:ibz1-ibz0+1
# 			pw[ibz0+iz-1,ibx1-ix+1] = boundary[2][iz,ix,it]
# 		end
# 	end
# end
# @fastmath @inbounds function boundaryf_t!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
# 	pw=p[:p]
# 	@simd for ix=1:ibx1-ibx0+1
# 		for iz=1:3
# 			pw[ibz0+iz-1,ibx0+ix-1] = boundary[1][iz,ix,it]
# 		end
# 	end
# end
# @fastmath @inbounds function boundaryf_b!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
# 	pw=p[:p]
# 	@simd for ix=1:ibx1-ibx0+1
# 		for iz=1:3
# 			pw[ibz1-iz+1,ibx0+ix-1] = boundary[3][iz,ix,it]
# 		end
# 	end
# end





# function get_boundary_indices(mgrid, ::FdtdOld)
#     # where to store the boundary values (careful, Born scatterers cannot be placed inside these!?
#     # while doing gradient tests)
#     # for 
#     ibx0 = _fd.npml - 2
#     ibx1 = length(mgrid[2]) + _fd.npml + 3
#     ibz0 = _fd.npml - 2
#     ibz1 = length(mgrid[1]) + _fd.npml + 3
#     #	ibx0=_fd.npml+1; ibx1=length(mgrid[2])+_fd.npml
#     #	ibz0=_fd.npml+1; ibz1=length(mgrid[1])+_fd.npml
#     #	println("**** Boundary Storage Changed **** ")

#     # for snaps
#     isx0, isz0 = _fd.npml, _fd.npml

#     return NamedArray(
#         [ibx0, ibx1, ibz0, ibz1, isx0, isz0],
#         ([:bx0, :bx1, :bz0, :bz1, :sx0, :sz0],),
#     )

# end

# function get_boundary_indices(mgrid, ::FdtdElastic)
#     return NamedArray([1], [:a]) # some dummy, update later	
# end

