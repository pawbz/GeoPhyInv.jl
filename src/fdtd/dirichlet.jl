

for dimnames in [zip([:1, :2, :3], [:z, :y, :x]), zip([:1, :2], [:z, :x])]
    is = broadcast(x -> Symbol(string("i", x)), getindex.(collect(dimnames), 2))
    ist = Meta.parse(string("(", [string(s, ",") for s in is]..., ")"))
    N = Meta.parse(string(length(is)))

    for (idim, dim) in dimnames
        i = Symbol("i", string(dim))
        # velocity-free boundary conditions    
        fname = Symbol("dirichlet", string(dim), "!")

        is1, is2, is3 = [replace(is, i => ii) for ii in [:1, :2, :3]]
        isn = replace(is, i => :n)
        isnm1 = replace(is, i => :(n - 1))
        isnp1 = replace(is, i => :(n + 1))

        irest = Meta.parse(
            string("(", [string(s, ",") for s in filter(x -> x != i, is)]..., ")"),
        )
        v = Symbol("v", string(dim))
        vs = broadcast(x -> Symbol(string("v", x)), getindex.(collect(dimnames), 2))
        vrest = filter(x -> x != v, vs)
        vrestt = Meta.parse(string("(", [string(s, ",") for s in vrest]..., ")"))

        @eval @parallel_indices($irest, function $fname($v::Data.Array{$N}, $(vrest...), n)
            # along other dimensions velocity grid matches tauii grid
            for vv in $vrestt
                vv[$(is1...)] = 0
                vv[$(isn...)] = 0
            end

            # relative to the tauii grid, velocity at i=2 should be zero, so making use of ghost cells 
            $v[$(is1...)] = -$v[$(is2...)]
            $v[$(isnp1...)] = -$v[$(isn...)]
            return
        end)
    end
end

