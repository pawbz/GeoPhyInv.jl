function lossvalue(m, loss, dobs::Records, pa::PFdtd, mparams=Medium(pa.c.attrib_mod))
    update!(pa, m, mparams)
    mode_save = pa.c.attrib_mod.mode
    pa.c.attrib_mod.mode = :forward
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
    pa.c.attrib_mod.mode = :forward
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
        end
    end
    # copy pac.gradients to g
    copyto!(g, Iterators.flatten(pa.c.gradients[mparams]))
    CUDA.allowscalar(false)
    pa.c.attrib_mod.mode = mode_save
    return lossvalue(loss, dobs, pa.c.data[1])
end