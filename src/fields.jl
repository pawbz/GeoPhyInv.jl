# Structs for stress and velocity fields (used for multiple dispatch)

abstract type Fields end

"""
Return a list of fields which contain `s`.
Retruns all the fields defined here by default.
"""
Fields(s = "") = filter(x -> contains(string(x), string(s)), nameof.(subtypes(Fields)))

# define pressure for acoustic 
struct p <: Fields end
# define stress and velocity (tauxx, tauxy, tauxz,)
for f in vcat(dim_names(3, "v"), dim_names(3, "tau"; order = 2))
    @eval struct $f <: Fields end
end
# define spatial derivatives of stress, velocity, pressure
dnames = dim_names(3, "d") # derivatives
vnames = dim_names(3, "v") # velocity
taunames = dim_names(3, "tau"; order = 2)
for (i, x) in enumerate(dim_names(3))
    tn = filter(y -> contains(string(y), string(x)), taunames)
    for t in vcat(tn, vnames, [:p])
        @eval struct $(Symbol("d", string(t), "d", string(x))) <: Fields end
    end
end

#-------------------------------------------------------
# grid indices and linear-interpolation weights
#-------------------------------------------------------
"""
Output CartesianIndices and corresponding linear interpolation weights for point P in mgrid
For example, when 2D,
mgrid=[range(1.3, stop=10.6,step=0.003), range(1.2,stop=15.3,step=0.004)]
P=[5,5]
"""
function get_indices_weights(mgrid, P)
    @assert length(mgrid) == length(P)
    N = length(mgrid)
    idx = [Interpolation.indminn(mgrid[i], P[i], 2) for i = 1:N]
    denomI = inv(prod([diff(mgrid[i][idx[i]])[1] for i = 1:N]))
    c = CartesianIndices(Tuple(broadcast(x -> x[1]:x[2], idx)))
    weights = zeros(length(c))
    for (i, cc) in enumerate(c)
        weights[i] = prod([abs(mgrid[i][cc[i]] - P[i]) for i = 1:N]) * denomI
    end
    return c, weights
end

for dimnames in [zip([:1, :2, :3], [:z, :y, :x]), zip([:1, :2], [:z, :x])]
    grids = broadcast(x -> Symbol(string("m", x)), getindex.(collect(dimnames), 2))
    for (idim, dim) in dimnames
        v = Symbol("v", string(dim))
        g = Symbol("m", string(dim))
        @eval function get_indices_weights(::$v, $(grids...), P) 
            # velocity grids are stagerred in respective dimensions, so generate half grids in respective dimensions 
            $g=range($g[1] - step($g) * (_fd.order - 1) * 0.5, step = step($g), length = length($g) + _fd.order - 1,)
            return get_indices_weights([$(grids...)], P)
        end
    end
end


#-------------------------------------------------------
# Acoustic
#-------------------------------------------------------
Base.zeros(::p, ::FdtdAcou, nz, nx) = Data.Array(zeros(nz, nx))
Base.zeros(::p, ::FdtdAcou, nz, ny, nx) = Data.Array(zeros(nz, ny, nx))
Base.zeros(::dpdx, ::FdtdAcou, nz, nx; pml::Bool = false) =
    Data.Array(zeros(nz - 2(_fd.order - 1), pml ? 2 * npml : nx - 1(_fd.order - 1)))
Base.zeros(::dpdx, ::FdtdAcou, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 2(_fd.order - 1),
        ny - 2(_fd.order - 1),
        pml ? 2 * npml : nx - 1(_fd.order - 1),
    ),
)
Base.zeros(::dpdz, ::FdtdAcou, nz, nx; pml::Bool = false) =
    Data.Array(zeros(pml ? 2 * npml : nz - 1(_fd.order - 1), nx - 2(_fd.order - 1)))
Base.zeros(::dpdz, ::FdtdAcou, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        pml ? 2 * npml : nz - 1(_fd.order - 1),
        ny - 2(_fd.order - 1),
        nx - 2(_fd.order - 1),
    ),
)
Base.zeros(::dpdy, ::FdtdAcou, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 2(_fd.order - 1),
        pml ? 2 * npml : ny - 1(_fd.order - 1),
        nx - 2(_fd.order - 1),
    ),
)
Base.zeros(::vx, ::FdtdAcou, nz, nx) = Data.Array(zeros(nz, nx + _fd.order - 1))
Base.zeros(::vz, ::FdtdAcou, nz, nx) = Data.Array(zeros(nz + _fd.order - 1, nx))
Base.zeros(::vx, ::FdtdAcou, nz, ny, nx) = Data.Array(zeros(nz, ny, nx + _fd.order - 1))
Base.zeros(::vz, ::FdtdAcou, nz, ny, nx) = Data.Array(zeros(nz + _fd.order - 1, ny, nx))
Base.zeros(::vy, ::FdtdAcou, nz, ny, nx) = Data.Array(zeros(nz, ny + 1, nx))
Base.zeros(::dvxdx, ::FdtdAcou, nz, nx; pml::Bool = false) =
    Data.Array(zeros(nz, pml ? 2 * npml : nx))
Base.zeros(::dvzdz, ::FdtdAcou, nz, nx; pml::Bool = false) =
    Data.Array(zeros(pml ? 2 * npml : nz, nx))
Base.zeros(::dvxdx, ::FdtdAcou, nz, ny, nx; pml::Bool = false) =
    Data.Array(zeros(nz, ny, pml ? 2 * npml : nx))
Base.zeros(::dvzdz, ::FdtdAcou, nz, ny, nx; pml::Bool = false) =
    Data.Array(zeros(pml ? 2 * npml : nz, ny, nx))
Base.zeros(::dvydy, ::FdtdAcou, nz, ny, nx; pml::Bool = false) =
    Data.Array(zeros(nz, pml ? 2 * npml : ny, nx))
Base.zeros(::dvxdz, ::FdtdAcou, nz, nx; pml::Bool = false) = Data.Array(zeros(1, 1))
Base.zeros(::dvzdx, ::FdtdAcou, nz, nx; pml::Bool = false) = Data.Array(zeros(1, 1))
Base.zeros(::dvxdz, ::FdtdAcou, nz, ny, nx; pml::Bool = false) = Data.Array(zeros(1, 1, 1))
Base.zeros(::dvydx, ::FdtdAcou, nz, ny, nx; pml::Bool = false) = Data.Array(zeros(1, 1, 1))
Base.zeros(::dvzdx, ::FdtdAcou, nz, ny, nx; pml::Bool = false) = Data.Array(zeros(1, 1, 1))
Base.zeros(::dvydz, ::FdtdAcou, nz, ny, nx; pml::Bool = false) = Data.Array(zeros(1, 1, 1))
Base.zeros(::dvxdy, ::FdtdAcou, nz, ny, nx; pml::Bool = false) = Data.Array(zeros(1, 1, 1))
Base.zeros(::dvzdy, ::FdtdAcou, nz, ny, nx; pml::Bool = false) = Data.Array(zeros(1, 1, 1))


#-------------------------------------------------------
# Elastic
#-------------------------------------------------------
Base.zeros(::tauxx, ::FdtdElastic, nz, ny, nx) = Data.Array(zeros(nz, ny, nx))
Base.zeros(::tauyy, ::FdtdElastic, nz, ny, nx) = Data.Array(zeros(nz, ny, nx))
Base.zeros(::tauzz, ::FdtdElastic, nz, ny, nx) = Data.Array(zeros(nz, ny, nx))
Base.zeros(::tauxy, ::FdtdElastic, nz, ny, nx) =
    Data.Array(zeros(nz - 2(_fd.order - 1), ny - 1(_fd.order - 1), nx - 1(_fd.order - 1)))
Base.zeros(::tauxz, ::FdtdElastic, nz, ny, nx) =
    Data.Array(zeros(nz - 1(_fd.order - 1), ny - 2(_fd.order - 1), nx - 1(_fd.order - 1)))
Base.zeros(::tauyz, ::FdtdElastic, nz, ny, nx) =
    Data.Array(zeros(nz - 1(_fd.order - 1), ny - 1(_fd.order - 1), nx - 2(_fd.order - 1)))

Base.zeros(::tauxx, ::FdtdElastic, nz, nx) = Data.Array(zeros(nz, nx))
Base.zeros(::tauzz, ::FdtdElastic, nz, nx) = Data.Array(zeros(nz, nx))
Base.zeros(::tauxz, ::FdtdElastic, nz, nx) =
    Data.Array(zeros(nz - 1(_fd.order - 1), nx - 1(_fd.order - 1)))

Base.zeros(::vx, ::FdtdElastic, nz, ny, nx) = Data.Array(zeros(nz, ny, nx + _fd.order - 1))
Base.zeros(::vy, ::FdtdElastic, nz, ny, nx) = Data.Array(zeros(nz, ny + _fd.order - 1, nx))
Base.zeros(::vz, ::FdtdElastic, nz, ny, nx) = Data.Array(zeros(nz + _fd.order - 1, ny, nx))

Base.zeros(::vx, ::FdtdElastic, nz, nx) = Data.Array(zeros(nz, nx + _fd.order - 1))
Base.zeros(::vz, ::FdtdElastic, nz, nx) = Data.Array(zeros(nz + _fd.order - 1, nx))

Base.zeros(::dvxdx, ::FdtdElastic, nz, ny, nx; pml::Bool = false) =
    Data.Array(zeros(nz, ny, pml ? 2 * npml : nx))
Base.zeros(::dvydy, ::FdtdElastic, nz, ny, nx; pml::Bool = false) =
    Data.Array(zeros(nz, pml ? 2 * npml : ny, nx))
Base.zeros(::dvzdz, ::FdtdElastic, nz, ny, nx; pml::Bool = false) =
    Data.Array(zeros(pml ? 2 * npml : nz, ny, nx))
Base.zeros(::dvxdy, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 2(_fd.order - 1),
        pml ? 2 * npml : ny - 1(_fd.order - 1),
        nx - 1(_fd.order - 1),
    ),
)
Base.zeros(::dvxdz, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        pml ? 2 * npml : nz - 1(_fd.order - 1),
        ny - 2(_fd.order - 1),
        nx - 1(_fd.order - 1),
    ),
)
Base.zeros(::dvydx, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 2(_fd.order - 1),
        ny - 1(_fd.order - 1),
        pml ? 2 * npml : nx - 1(_fd.order - 1),
    ),
)
Base.zeros(::dvydz, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        pml ? 2 * npml : nz - 1(_fd.order - 1),
        ny - 1(_fd.order - 1),
        nx - 2(_fd.order - 1),
    ),
)
Base.zeros(::dvzdx, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 1(_fd.order - 1),
        ny - 2(_fd.order - 1),
        pml ? 2 * npml : nx - 1(_fd.order - 1),
    ),
)
Base.zeros(::dvzdy, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 1(_fd.order - 1),
        pml ? 2 * npml : ny - 1(_fd.order - 1),
        nx - 2(_fd.order - 1),
    ),
)

Base.zeros(::dvxdx, ::FdtdElastic, nz, nx; pml::Bool = false) =
    Data.Array(zeros(nz, pml ? 2 * npml : nx))
Base.zeros(::dvzdz, ::FdtdElastic, nz, nx; pml::Bool = false) =
    Data.Array(zeros(pml ? 2 * npml : nz, nx))
Base.zeros(::dvxdz, ::FdtdElastic, nz, nx; pml::Bool = false) =
    Data.Array(zeros(pml ? 2 * npml : nz - 1(_fd.order - 1), nx - 1(_fd.order - 1)))
Base.zeros(::dvzdx, ::FdtdElastic, nz, nx; pml::Bool = false) =
    Data.Array(zeros(nz - 1(_fd.order - 1), pml ? 2 * npml : nx - 1(_fd.order - 1)))

Base.zeros(::dtauxxdx, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 2(_fd.order - 1),
        ny - 2(_fd.order - 1),
        pml ? 2 * npml : nx - 1(_fd.order - 1),
    ),
)
Base.zeros(::dtauyydy, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 2(_fd.order - 1),
        pml ? 2 * npml : ny - 1(_fd.order - 1),
        nx - 2(_fd.order - 1),
    ),
)
Base.zeros(::dtauzzdz, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        pml ? 2 * npml : nz - 1(_fd.order - 1),
        ny - 2(_fd.order - 1),
        nx - 2(_fd.order - 1),
    ),
)
Base.zeros(::dtauxydx, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 2(_fd.order - 1),
        ny - 1(_fd.order - 1),
        pml ? 2 * npml : nx - 2(_fd.order - 1),
    ),
)
Base.zeros(::dtauxydy, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 2(_fd.order - 1),
        pml ? 2 * npml : ny - 2(_fd.order - 1),
        nx - 1(_fd.order - 1),
    ),
)
Base.zeros(::dtauxzdx, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 1(_fd.order - 1),
        ny - 2(_fd.order - 1),
        pml ? 2 * npml : nx - 2(_fd.order - 1),
    ),
)
Base.zeros(::dtauxzdz, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        pml ? 2 * npml : nz - 2(_fd.order - 1),
        ny - 2(_fd.order - 1),
        nx - 1(_fd.order - 1),
    ),
)
Base.zeros(::dtauyzdy, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 1(_fd.order - 1),
        pml ? 2 * npml : ny - 2(_fd.order - 1),
        nx - 2(_fd.order - 1),
    ),
)
Base.zeros(::dtauyzdz, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        pml ? 2 * npml : nz - 2(_fd.order - 1),
        ny - 1(_fd.order - 1),
        nx - 2(_fd.order - 1),
    ),
)

Base.zeros(::dtauxxdx, ::FdtdElastic, nz, nx; pml::Bool = false) =
    Data.Array(zeros(nz - 2(_fd.order - 1), pml ? 2 * npml : nx - 1(_fd.order - 1)))
Base.zeros(::dtauzzdz, ::FdtdElastic, nz, nx; pml::Bool = false) =
    Data.Array(zeros(pml ? 2 * npml : nz - 1(_fd.order - 1), nx - 2(_fd.order - 1)))
Base.zeros(::dtauxzdx, ::FdtdElastic, nz, nx; pml::Bool = false) =
    Data.Array(zeros(nz - 1(_fd.order - 1), pml ? 2 * npml : nx - 2(_fd.order - 1)))
Base.zeros(::dtauxzdz, ::FdtdElastic, nz, nx; pml::Bool = false) =
    Data.Array(zeros(pml ? 2 * npml : nz - 2(_fd.order - 1), nx - 1(_fd.order - 1)))

# Dont need pressure for FdtdElastic, so no initialization
for f in Fields("p")
    @eval Base.zeros(::$f, ::FdtdElastic, args1...; args2...) =
        Data.Array(zeros(fill(1, length(args1))...))
end
# remove stress from acoustic simulations
for f in Fields("tau")
    @eval Base.zeros(::$f, ::FdtdAcou, args1...; args2...) =
        Data.Array(zeros(fill(1, length(args1))...))
end
# no "y" for 2D simulations
for f in Fields("y")
    @eval Base.zeros(::$f, attrib_mod, nz, nx; args...) = Data.Array(zeros(fill(1, 2)...))
end


