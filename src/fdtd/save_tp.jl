
"""
Save the wavefield to `tp` before moving to the next time step (used for gradient calculation only)
"""
function save_tp!(::Val{:adjoint}, pap)
    w1t = pap.w1[:t]
    w1tp = pap.w1[:tp]
    broadcast(names(w1tp)[1]) do field
        copyto!(w1tp[field], w1t[field])
    end
end


function save_tp!(::Any, pap)
end