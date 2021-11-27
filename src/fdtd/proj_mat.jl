
"""
Output CartesianIndices and corresponding linear interpolation weights for point P in mgrid
For example, when 2D,
mgrid=[range(1.3, stop=10.6,step=0.003), range(1.2,stop=15.3,step=0.004)]
P=[5,5]
"""
function update_proj_mat!(R, pac, sr, f, mgrid, P)
    @assert length(mgrid) == length(P)
    N = length(mgrid)
    idx = [Interpolation.indminn(mgrid[i], P[i], 2) for i = 1:N]
    denomI = inv(prod([diff(mgrid[i][idx[i]])[1] for i = 1:N]))
    c = CartesianIndices(Tuple(broadcast(x -> x[1]:x[2], idx)))
    @views RR = reshape(R, length.(mgrid)...)
    for (i, cc) in enumerate(c)
        RR[cc] =
            weight(sr, f, prod([abs(mgrid[i][cc[i]] - P[i]) for i = 1:N]) * denomI, pac)
    end
    return R
end

for dimnames in [zip([:1, :2, :3], [:z, :y, :x]), zip([:1, :2], [:z, :x])]
    grids = broadcast(x -> Symbol(string("m", x)), getindex.(collect(dimnames), 2))
    # get indices for any fields, routed via mgrid adjustment
    @eval function update_proj_mat!(R, pac, sr, f::Any, attrib_mod, $(grids...), P)
        return update_proj_mat!(R, pac, sr, f, get_mgrid(f, attrib_mod, $(grids...)), P)
    end
end


# multiplication with modK
# division of source term with δx and δz (see Jan's fdelmodc manual)
# on pressure grid
# * pac.mod[:K][si] * pac.fc[:dt] * prod(inv.(step.(pac.medium.mgrid)))
weight(::Srcs, ::Any, w, pac) =
    w * pac.fc[:dt] * prod(inv.(step.(pac.medium.mgrid)))

weight(::Recs, ::Any, w, pac) = w
