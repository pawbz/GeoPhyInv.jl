function lossvalue(m, loss, dobs::Records, pa::PFdtd)
    update!(pa, m)
    mode_save = pa.c.attrib_mod.mode
    pa.c.attrib_mod.mode = :forward
    update!(pa, pa.c.srcwav, [1, 0], verbose=false)
    update!(pa)
    pa.c.attrib_mod.mode = mode_save
    return lossvalue(loss, dobs, pa.c.data[1])
end

function gradient!(g, m, loss, dobs::Records, pa::PFdtd)
    # put m into pac.mod
    update!(pa, m)

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
    mchunks = chunk(m, size=length(first(pa.c.mod)))
    foreach(pa.c.gradients, mchunks, pa.c.ref_mod) do gm, m, r
        map!(gm, gm, m) do gm1, m1
            # chainrule
            gm1 * exp(m1) * r
        end
    end
    # copy pac.gradients to g
    copyto!(g, Iterators.flatten(pa.c.gradients))
    CUDA.allowscalar(false)
    pa.c.attrib_mod.mode = mode_save
    return lossvalue(loss, dobs, pa.c.data[1])
end