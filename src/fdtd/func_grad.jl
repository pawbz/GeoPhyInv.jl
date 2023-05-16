function lossvalue(m, loss, dobs::Records, pa::PFdtd, mparams=Medium(pa.c.attrib_mod))
    update!(pa, m, mparams)
    mode_save = pa.c.attrib_mod.mode
    pa.c.attrib_mod.mode = :forward_save
    update!(pa, pa.c.srcwav, [1, 0], verbose=false)
    update!(pa)
    pa.c.attrib_mod.mode = mode_save
    return lossvalue(loss, dobs, pa.c.data[1])
end

function gradient!(g, m, loss, dobs::Records, pa::PFdtd, mparams=Medium(pa.c.attrib_mod))

    # put m into pac.mod
    update!(pa, m, mparams)

    mode_save = pa.c.attrib_mod.mode

    # generate data and save boundaries
    pa.c.attrib_mod.mode = :forward_save
    update!(pa, pa.c.srcwav, [1, 0], verbose=false)
    update!(pa)

    # adjoint sources
    gradient!(pa.c.srcwav[2], loss, dobs, pa.c.data[1])
    reverse!(pa.c.srcwav[2])
    update!(pa, pa.c.srcwav, [-1, 1], verbose=false)
    pa.c.attrib_mod.mode = :adjoint
    update!(pa)

    CUDA.allowscalar(true)
    mchunks = chunk(m, length(mparams))
    broadcast(enumerate(mparams)) do (i, mname)
        r = pa.c.ref_mod[mname]
        gm = pa.c.gradients[mname]
        x = mchunks[i]
        # chainrule
        map!(gm, gm, x) do gm1, x1
            gm1 * exp(x1) * r
            # gm1 * r
        end
    end
    # copy pac.gradients to g
    copyto!(g, Iterators.flatten(pa.c.gradients[mparams]))
    CUDA.allowscalar(false)
    pa.c.attrib_mod.mode = mode_save
    return lossvalue(loss, dobs, pa.c.data[1])
end

# ============== Born ========================================
# m is nondimensionalized model vector
function forward_map!(d, m, pa::PFdtd)
    # copy input m to pac.δmod
    broadcast(pa.c.δmod, chunk(m, size=length(first(pa.c.δmod)))) do m1, m2
        CUDA.@allowscalar copyto!(m1, m2)
    end

    update!(pa, pa.c.srcwav, [1, 1], verbose=false)
    
    mode_save = pa.c.attrib_mod.mode
    pa.c.attrib_mod.mode = :forward
    update!(pa)
    # copy pac.data to d (only implemented for first supersource for now
    copyto!(d, Iterators.flatten(pa.c.data[2][1].d))
    pa.c.attrib_mod.mode = mode_save
end

function adjoint_map!(gm, d, pa::PFdtd)

    # copy input d to pac.srcwav (only implemented first source for now
    broadcast(pa.c.srcwav[2][1].d, chunk(d, size=length(first(pa.c.srcwav[2][1].d)))) do d1, d2
        CUDA.@allowscalar copyto!(d1, d2)
    end

    # time reversal
    foreach(pa.c.srcwav[2]) do S # each supersource
        foreach(S.d) do s # each field
            reverse!(s, dims=1)
        end
    end

    update!(pa, pa.c.srcwav, [-1, 1], verbose=false)

    mode_save = pa.c.attrib_mod.mode
    pa.c.attrib_mod.mode = :adjoint
    update!(pa)

    # copy pac.gradients to m
    copyto!(gm, Iterators.flatten(pa.c.gradients))
    pa.c.attrib_mod.mode = mode_save
end



"""
```
F=LinearMap(pa)
```
If `pa` is an instance of `SeisInvExpt`, then 
return the linearized forward modeling operator `F`, such that
`F*x` can be computed without explicitly storing the operator matrix (see `LinearMaps.jl`).
The imaging/migration operator is given by `transpose(F)`. 
These operators are the building blocks of iterative optimization schemes.
"""
function LinearMaps.LinearMap(pa::T) where {T<:PFdtd{FdtdAcoustic{Born}}}
    fw = (y, x) -> forward_map!(y, x, pa)
    bk = (y, x) -> adjoint_map!(y, x, pa)

    # data (output) length
    nd = mapreduce(+, pa.c.data[1]) do d
        return sum(length.(d.d))
    end
    # length of medium (input)
    nm = sum(length.(pa.c.δmod))

    return LinearMap(fw, bk,
        nd, nm, 
        ismutating=true)
end


