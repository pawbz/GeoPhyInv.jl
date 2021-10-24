

for dimnames in [zip([:1, :2, :3], [:z, :y, :x]), zip([:1, :2], [:z, :x])]
    is = broadcast(x -> Symbol(string("i", x)), getindex.(collect(dimnames), 2))
    ist = Meta.parse(string("(", [string(s, ",") for s in is]..., ")"))
    N = Meta.parse(string(length(is)))

    for (idim, dim) in dimnames
        i = Symbol("i", string(dim))
        # velocity-free boundary conditions    
        fname = Symbol("dirichlet", string(dim), "!")

        fdh=div(_fd.order,2)
        # mirror ghost cells 
        ighost=vcat(
            [[replace(is, i => :($ifd)), replace(is, i => :($(_fd.order+1-ifd)))] for ifd = 1:fdh],
            [[replace(is, i => :(n + $(_fd.order - ifd))), replace(is, i => :(n + $(ifd-1)))] for ifd = 1:fdh]
        )

        isn = replace(is, i => :n)
        is1 = replace(is, i => :1)

        irest = Meta.parse(
            string("(", [string(s, ",") for s in filter(x -> x != i, is)]..., ")"),
        )
        v = Symbol("v", string(dim))
        vs = broadcast(x -> Symbol(string("v", x)), getindex.(collect(dimnames), 2))
        vrest = filter(x -> x != v, vs)

        @eval @parallel_indices($irest, function $fname($v::Data.Array{$N}, $(vrest...), n)
            # along other dimensions velocity grid matches tauii grid
            $((quote $vv[$(is1...)] = 0; $vv[$(isn...)] = 0 end for vv in vrest)...)
            # relative to the tauii grid, velocity at i=2 should be zero, so making use of ghost cells 
            $((quote $v[$(ig[1]...)] = -$v[$(ig[2]...)] end for ig in ighost)...)
            return
        end)
    end
end

