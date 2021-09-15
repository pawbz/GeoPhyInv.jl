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
Fields(s = "") = filter(x -> contains(string(x), string(s)), subtypes(Fields))

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


Base.zeros(::tauxx, ::FdtdElastic, nz, ny, nx) = Data.Array(zeros(nz, ny, nx))
Base.zeros(::tauyy, ::FdtdElastic, nz, ny, nx) = Data.Array(zeros(nz, ny, nx))
Base.zeros(::tauzz, ::FdtdElastic, nz, ny, nx) = Data.Array(zeros(nz, ny, nx))
Base.zeros(::tauxy, ::FdtdElastic, nz, ny, nx) = Data.Array(zeros(nz - 2, ny - 1, nx - 1))
Base.zeros(::tauxz, ::FdtdElastic, nz, ny, nx) = Data.Array(zeros(nz - 1, ny - 2, nx - 1))
Base.zeros(::tauyz, ::FdtdElastic, nz, ny, nx) = Data.Array(zeros(nz - 1, ny - 1, nx - 2))

Base.zeros(::vx, ::FdtdElastic, nz, ny, nx) = Data.Array(zeros(nz, ny, nx + 1))
Base.zeros(::vy, ::FdtdElastic, nz, ny, nx) = Data.Array(zeros(nz, ny + 1, nx))
Base.zeros(::vz, ::FdtdElastic, nz, ny, nx) = Data.Array(zeros(nz + 1, ny, nx))

Base.zeros(::dvxdx, ::FdtdElastic, nz, ny, nx; pml::Bool=false) = Data.Array(zeros(nz, ny, pml ? 2*npml : nx))
Base.zeros(::dvydy, ::FdtdElastic, nz, ny, nx; pml::Bool=false) = Data.Array(zeros(nz, pml ? 2*npml : ny, nx))
Base.zeros(::dvzdz, ::FdtdElastic, nz, ny, nx; pml::Bool=false) = Data.Array(zeros(pml ? 2*npml : nz, ny, nx))
Base.zeros(::dvxdy, ::FdtdElastic, nz, ny, nx; pml::Bool=false) = Data.Array(zeros(nz - 2, pml ? 2*npml : ny - 1, nx - 1))
Base.zeros(::dvxdz, ::FdtdElastic, nz, ny, nx; pml::Bool=false) = Data.Array(zeros(pml ? 2*npml : nz - 1, ny - 2, nx - 1))
Base.zeros(::dvydx, ::FdtdElastic, nz, ny, nx; pml::Bool=false) = Data.Array(zeros(nz - 2, ny - 1, pml ? 2*npml : nx - 1))
Base.zeros(::dvydz, ::FdtdElastic, nz, ny, nx; pml::Bool=false) = Data.Array(zeros(pml ? 2*npml : nz - 1, ny - 1, nx - 2))
Base.zeros(::dvzdx, ::FdtdElastic, nz, ny, nx; pml::Bool=false) = Data.Array(zeros(nz - 1, ny - 2, pml ? 2*npml : nx - 1))
Base.zeros(::dvzdy, ::FdtdElastic, nz, ny, nx; pml::Bool=false) = Data.Array(zeros(nz - 1, pml ? 2*npml : ny - 1, nx - 2))

Base.zeros(::dtauxxdx, ::FdtdElastic, nz, ny, nx; pml::Bool=false) =
    Data.Array(zeros(nz - 2, ny - 2, pml ? 2*npml : nx - 1))
Base.zeros(::dtauyydy, ::FdtdElastic, nz, ny, nx; pml::Bool=false) =
    Data.Array(zeros(nz - 2, pml ? 2*npml : ny - 1, nx - 2))
Base.zeros(::dtauzzdz, ::FdtdElastic, nz, ny, nx; pml::Bool=false) =
    Data.Array(zeros(pml ? 2*npml : nz - 1, ny - 2, nx - 2))
Base.zeros(::dtauxydx, ::FdtdElastic, nz, ny, nx; pml::Bool=false) =
    Data.Array(zeros(nz - 2, ny - 1, pml ? 2*npml : nx - 2))
Base.zeros(::dtauxydy, ::FdtdElastic, nz, ny, nx; pml::Bool=false) =
    Data.Array(zeros(nz - 2, pml ? 2*npml : ny - 2, nx - 1))
Base.zeros(::dtauxzdx, ::FdtdElastic, nz, ny, nx; pml::Bool=false) =
    Data.Array(zeros(nz - 1, ny - 2, pml ? 2*npml : nx - 2))
Base.zeros(::dtauxzdz, ::FdtdElastic, nz, ny, nx; pml::Bool=false) =
    Data.Array(zeros(pml ? 2*npml : nz - 2, ny - 2, nx - 1))
Base.zeros(::dtauyzdy, ::FdtdElastic, nz, ny, nx; pml::Bool=false) =
    Data.Array(zeros(nz - 1, pml ? 2*npml : ny - 2, nx - 2))
Base.zeros(::dtauyzdz, ::FdtdElastic, nz, ny, nx; pml::Bool=false) =
    Data.Array(zeros(pml ? 2*npml : nz - 2, ny - 1, nx - 2))

# Dont need pressure for FdtdElastic, so no initialization
for f in Fields("p")
    @eval Base.zeros(::$f, ::FdtdElastic, args1...; args2...) = Data.Array(zeros(fill(1, length(args1))...))
end
# remove stress 
# for f in Fields("tau")
    # Base.zeros(::f, ::FdtdAcou, args...) = Data.Array(zeros(fill(1, length(args))...))
# end


