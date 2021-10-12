# Structs for stress and velocity fields (used for multiple dispatch)

"""
Return axis names of 1D, 2D or 3D fields
"""
function dim_names(N, prefix = "", suffix = ""; order = 1)
    if (N == 3)
        names = (order == 2) ? [:zz, :yy, :xx, :xz, :xy, :yz] : [:z, :y, :x]
    elseif (N == 2)
        names = (order == 2) ? [:zz, :xx, :xz] : [:z, :x]
    elseif (N == 1)
        names = [:z]
    else
        error("invalid dim num")
    end
    return [Symbol(prefix, string(m), suffix) for m in names]
end

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
# Acoustic
#-------------------------------------------------------
Base.zeros(::p, ::FdtdAcou, nz, nx) = Data.Array(zeros(nz, nx))
Base.zeros(::p, ::FdtdAcou, nz, ny, nx) = Data.Array(zeros(nz, ny, nx))
Base.zeros(::dpdx, ::FdtdAcou, nz, nx; pml::Bool = false) =
    Data.Array(zeros(nz - 2(FD_ORDER - 1), pml ? 2 * npml : nx - 1(FD_ORDER - 1)))
Base.zeros(::dpdx, ::FdtdAcou, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 2(FD_ORDER - 1),
        ny - 2(FD_ORDER - 1),
        pml ? 2 * npml : nx - 1(FD_ORDER - 1),
    ),
)
Base.zeros(::dpdz, ::FdtdAcou, nz, nx; pml::Bool = false) =
    Data.Array(zeros(pml ? 2 * npml : nz - 1(FD_ORDER - 1), nx - 2(FD_ORDER - 1)))
Base.zeros(::dpdz, ::FdtdAcou, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        pml ? 2 * npml : nz - 1(FD_ORDER - 1),
        ny - 2(FD_ORDER - 1),
        nx - 2(FD_ORDER - 1),
    ),
)
Base.zeros(::dpdy, ::FdtdAcou, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 2(FD_ORDER - 1),
        pml ? 2 * npml : ny - 1(FD_ORDER - 1),
        nx - 2(FD_ORDER - 1),
    ),
)
Base.zeros(::vx, ::FdtdAcou, nz, nx) = Data.Array(zeros(nz, nx + 1))
Base.zeros(::vz, ::FdtdAcou, nz, nx) = Data.Array(zeros(nz + 1, nx))
Base.zeros(::vx, ::FdtdAcou, nz, ny, nx) = Data.Array(zeros(nz, ny, nx + 1))
Base.zeros(::vz, ::FdtdAcou, nz, ny, nx) = Data.Array(zeros(nz + 1, ny, nx))
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
    Data.Array(zeros(nz - 2(FD_ORDER - 1), ny - 1(FD_ORDER - 1), nx - 1(FD_ORDER - 1)))
Base.zeros(::tauxz, ::FdtdElastic, nz, ny, nx) =
    Data.Array(zeros(nz - 1(FD_ORDER - 1), ny - 2(FD_ORDER - 1), nx - 1(FD_ORDER - 1)))
Base.zeros(::tauyz, ::FdtdElastic, nz, ny, nx) =
    Data.Array(zeros(nz - 1(FD_ORDER - 1), ny - 1(FD_ORDER - 1), nx - 2(FD_ORDER - 1)))

Base.zeros(::tauxx, ::FdtdElastic, nz, nx) = Data.Array(zeros(nz, nx))
Base.zeros(::tauzz, ::FdtdElastic, nz, nx) = Data.Array(zeros(nz, nx))
Base.zeros(::tauxz, ::FdtdElastic, nz, nx) =
    Data.Array(zeros(nz - 1(FD_ORDER - 1), nx - 1(FD_ORDER - 1)))

Base.zeros(::vx, ::FdtdElastic, nz, ny, nx) = Data.Array(zeros(nz, ny, nx + FD_ORDER - 1))
Base.zeros(::vy, ::FdtdElastic, nz, ny, nx) = Data.Array(zeros(nz, ny + FD_ORDER - 1, nx))
Base.zeros(::vz, ::FdtdElastic, nz, ny, nx) = Data.Array(zeros(nz + FD_ORDER - 1, ny, nx))

Base.zeros(::vx, ::FdtdElastic, nz, nx) = Data.Array(zeros(nz, nx + 1))
Base.zeros(::vz, ::FdtdElastic, nz, nx) = Data.Array(zeros(nz + 1, nx))

Base.zeros(::dvxdx, ::FdtdElastic, nz, ny, nx; pml::Bool = false) =
    Data.Array(zeros(nz, ny, pml ? 2 * npml : nx))
Base.zeros(::dvydy, ::FdtdElastic, nz, ny, nx; pml::Bool = false) =
    Data.Array(zeros(nz, pml ? 2 * npml : ny, nx))
Base.zeros(::dvzdz, ::FdtdElastic, nz, ny, nx; pml::Bool = false) =
    Data.Array(zeros(pml ? 2 * npml : nz, ny, nx))
Base.zeros(::dvxdy, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 2(FD_ORDER - 1),
        pml ? 2 * npml : ny - 1(FD_ORDER - 1),
        nx - 1(FD_ORDER - 1),
    ),
)
Base.zeros(::dvxdz, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        pml ? 2 * npml : nz - 1(FD_ORDER - 1),
        ny - 2(FD_ORDER - 1),
        nx - 1(FD_ORDER - 1),
    ),
)
Base.zeros(::dvydx, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 2(FD_ORDER - 1),
        ny - 1(FD_ORDER - 1),
        pml ? 2 * npml : nx - 1(FD_ORDER - 1),
    ),
)
Base.zeros(::dvydz, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        pml ? 2 * npml : nz - 1(FD_ORDER - 1),
        ny - 1(FD_ORDER - 1),
        nx - 2(FD_ORDER - 1),
    ),
)
Base.zeros(::dvzdx, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 1(FD_ORDER - 1),
        ny - 2(FD_ORDER - 1),
        pml ? 2 * npml : nx - 1(FD_ORDER - 1),
    ),
)
Base.zeros(::dvzdy, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 1(FD_ORDER - 1),
        pml ? 2 * npml : ny - 1(FD_ORDER - 1),
        nx - 2(FD_ORDER - 1),
    ),
)

Base.zeros(::dvxdx, ::FdtdElastic, nz, nx; pml::Bool = false) =
    Data.Array(zeros(nz, pml ? 2 * npml : nx))
Base.zeros(::dvzdz, ::FdtdElastic, nz, nx; pml::Bool = false) =
    Data.Array(zeros(pml ? 2 * npml : nz, nx))
Base.zeros(::dvxdz, ::FdtdElastic, nz, nx; pml::Bool = false) =
    Data.Array(zeros(pml ? 2 * npml : nz - 1(FD_ORDER - 1), nx - 1(FD_ORDER - 1)))
Base.zeros(::dvzdx, ::FdtdElastic, nz, nx; pml::Bool = false) =
    Data.Array(zeros(nz - 1(FD_ORDER - 1), pml ? 2 * npml : nx - 1(FD_ORDER - 1)))

Base.zeros(::dtauxxdx, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 2(FD_ORDER - 1),
        ny - 2(FD_ORDER - 1),
        pml ? 2 * npml : nx - 1(FD_ORDER - 1),
    ),
)
Base.zeros(::dtauyydy, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 2(FD_ORDER - 1),
        pml ? 2 * npml : ny - 1(FD_ORDER - 1),
        nx - 2(FD_ORDER - 1),
    ),
)
Base.zeros(::dtauzzdz, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        pml ? 2 * npml : nz - 1(FD_ORDER - 1),
        ny - 2(FD_ORDER - 1),
        nx - 2(FD_ORDER - 1),
    ),
)
Base.zeros(::dtauxydx, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 2(FD_ORDER - 1),
        ny - 1(FD_ORDER - 1),
        pml ? 2 * npml : nx - 2(FD_ORDER - 1),
    ),
)
Base.zeros(::dtauxydy, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 2(FD_ORDER - 1),
        pml ? 2 * npml : ny - 2(FD_ORDER - 1),
        nx - 1(FD_ORDER - 1),
    ),
)
Base.zeros(::dtauxzdx, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 1(FD_ORDER - 1),
        ny - 2(FD_ORDER - 1),
        pml ? 2 * npml : nx - 2(FD_ORDER - 1),
    ),
)
Base.zeros(::dtauxzdz, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        pml ? 2 * npml : nz - 2(FD_ORDER - 1),
        ny - 2(FD_ORDER - 1),
        nx - 1(FD_ORDER - 1),
    ),
)
Base.zeros(::dtauyzdy, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        nz - 1(FD_ORDER - 1),
        pml ? 2 * npml : ny - 2(FD_ORDER - 1),
        nx - 2(FD_ORDER - 1),
    ),
)
Base.zeros(::dtauyzdz, ::FdtdElastic, nz, ny, nx; pml::Bool = false) = Data.Array(
    zeros(
        pml ? 2 * npml : nz - 2(FD_ORDER - 1),
        ny - 1(FD_ORDER - 1),
        nx - 2(FD_ORDER - 1),
    ),
)

Base.zeros(::dtauxxdx, ::FdtdElastic, nz, nx; pml::Bool = false) =
    Data.Array(zeros(nz - 2(FD_ORDER - 1), pml ? 2 * npml : nx - 1(FD_ORDER - 1)))
Base.zeros(::dtauzzdz, ::FdtdElastic, nz, nx; pml::Bool = false) =
    Data.Array(zeros(pml ? 2 * npml : nz - 1(FD_ORDER - 1), nx - 2(FD_ORDER - 1)))
Base.zeros(::dtauxzdx, ::FdtdElastic, nz, nx; pml::Bool = false) =
    Data.Array(zeros(nz - 1(FD_ORDER - 1), pml ? 2 * npml : nx - 2(FD_ORDER - 1)))
Base.zeros(::dtauxzdz, ::FdtdElastic, nz, nx; pml::Bool = false) =
    Data.Array(zeros(pml ? 2 * npml : nz - 2(FD_ORDER - 1), nx - 1(FD_ORDER - 1)))

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


