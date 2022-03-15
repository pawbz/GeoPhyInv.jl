# Structs for stress and velocity fields (used for multiple dispatch)

abstract type Fields end

"""
Return a list of fields which contain `s`.
Returns all the elastic fields defined here by default i.e., for `ndims=3`.
"""
function Fields(s = ""; ndims = _fd.ndims)
    @assert(ndims ∈ [2, 3])
    f = filter(x -> contains(string(x), string(s)), nameof.(subtypes(Fields)))
    if (ndims == 2)
        f = filter(x -> !contains(string(x), "y"), f) # filter out y when 2D
    end
    return f
end
function Fields(::T, s = ""; ndims = _fd.ndims) where {T<:Elastic}
    f = Fields(s, ndims = ndims) # get all
    f = filter(x -> !contains(string(x), "p"), f) # filter out pressure
end
function Fields(::T, s = ""; ndims = _fd.ndims) where {T<:Acoustic}
    f = Fields(s, ndims = ndims) # get all
    f = filter(x -> !contains(string(x), "tau"), f) # filter out tau
    f = filter(x -> !in(x, [:dvxdy, :dvxdz, :dvydx, :dvydz, :dvzdx, :dvzdy]), f) # filter out dvxdy
end

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

for dimnames in [zip([:1, :2, :3], [:z, :y, :x]), zip([:1, :2], [:z, :x])]
    grids = broadcast(x -> Symbol(string("m", x)), getindex.(collect(dimnames), 2))
    sizes1 = broadcast(x -> Symbol(string("n1", x)), getindex.(collect(dimnames), 2))
    sizes = broadcast(x -> Symbol(string("n", x)), getindex.(collect(dimnames), 2))
    # initialize fields by simply getting grids 
    for f in Fields(ndims = length(collect(dimnames)))
        n1pml = Symbol("n1", string(last(string(f))))
        @eval function Base.zeros(
            ::$f,
            attrib_mod,
            $(sizes...);
            pml::Bool = false,
        ) where {T<:Fields}
            # get new sizes
            $(
                (
                    quote
                        $n1 = length(
                            get_mgrid(
                                $f(),
                                attrib_mod,
                                $((
                                        quote
                                            range(0, length = $n, stop = 1)
                                        end for n in sizes
                                    )...),
                            )[$i],
                        )
                    end for (i, n1) in enumerate(sizes1)
                )...
            )
            # if pml flag size on corresponding dimension is just 2*npml
            if (pml)
                $n1pml = 2 * _fd.npml
            end
            return Data.Array(zeros($(sizes1...)))
        end
    end
end


#-------------------------------------------------------
# Acoustic
#-------------------------------------------------------
get_mgrid(::p, ::FdtdAcoustic, mz, mx) = [mz, mx]
get_mgrid(::p, ::FdtdAcoustic, mz, my, mx) = [mz, my, mx]
get_mgrid(::dpdx, ::FdtdAcoustic, mz, mx) = [
    range(
        mz[1] + step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - 2(_fd.order - 1),
    ),
    range(
        mx[1] + 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - (_fd.order - 1),
    ),
]
get_mgrid(::dpdx, ::FdtdAcoustic, mz, my, mx) = [
    range(
        mz[1] + step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - 2(_fd.order - 1),
    ),
    range(
        my[1] + step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - 2(_fd.order - 1),
    ),
    range(
        mx[1] + 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - (_fd.order - 1),
    ),
]
get_mgrid(::dpdz, ::FdtdAcoustic, mz, mx) = [
    range(
        mz[1] + 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - (_fd.order - 1),
    ),
    range(
        mx[1] + step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - 2(_fd.order - 1),
    ),
]
get_mgrid(::dpdz, ::FdtdAcoustic, mz, my, mx) = [
    range(
        mz[1] + 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - (_fd.order - 1),
    ),
    range(
        my[1] + step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - 2(_fd.order - 1),
    ),
    range(
        mx[1] + step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - 2(_fd.order - 1),
    ),
]
get_mgrid(::dpdy, ::FdtdAcoustic, mz, my, mx) = [
    range(
        mz[1] + step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - 2(_fd.order - 1),
    ),
    range(
        my[1] + 0.5 * step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - (_fd.order - 1),
    ),
    range(
        mx[1] + step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - 2(_fd.order - 1),
    ),
]
get_mgrid(::vx, ::FdtdAcoustic, mz, mx) = [
    mz,
    range(
        mx[1] - 0.5 *  step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) + (_fd.order - 1),
    ),
]
get_mgrid(::vz, ::FdtdAcoustic, mz, mx) = [
    range(
        mz[1] - 0.5 * step(mz) *(_fd.order - 1),
        step = step(mz),
        length = length(mz) + (_fd.order - 1),
    ),
    mx,
]
get_mgrid(::vx, ::FdtdAcoustic, mz, my, mx) = [
    mz,
    my,
    range(
        mx[1] - 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) + (_fd.order - 1),
    ),
]
get_mgrid(::vz, ::FdtdAcoustic, mz, my, mx) = [
    range(
        mz[1] - 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) + (_fd.order - 1),
    ),
    my,
    mx,
]
get_mgrid(::vy, ::FdtdAcoustic, mz, my, mx) = [
    mz,
    range(
        my[1] - 0.5 * step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) + (_fd.order - 1),
    ),
    mx,
]
get_mgrid(::dvxdx, ::FdtdAcoustic, mz, mx) = [mz, mx]
get_mgrid(::dvzdz, ::FdtdAcoustic, mz, mx) = [mz, mx]
get_mgrid(::dvxdx, ::FdtdAcoustic, mz, my, mx) = [mz, my, mx]
get_mgrid(::dvzdz, ::FdtdAcoustic, mz, my, mx) = [mz, my, mx]
get_mgrid(::dvydy, ::FdtdAcoustic, mz, my, mx) = [mz, my, mx]


#-------------------------------------------------------
# Elastic
#-------------------------------------------------------
get_mgrid(::tauxx, ::FdtdElastic, mz, my, mx) = [mz, my, mx]
get_mgrid(::tauyy, ::FdtdElastic, mz, my, mx) = [mz, my, mx]
get_mgrid(::tauzz, ::FdtdElastic, mz, my, mx) = [mz, my, mx]
get_mgrid(::tauxy, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] + step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - 2(_fd.order - 1),
    ),
    range(
        my[1] + 0.5 * step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - (_fd.order - 1),
    ),
    range(
        mx[1] + 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - (_fd.order - 1),
    ),
]
get_mgrid(::tauxz, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] + 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - (_fd.order - 1),
    ),
    range(
        my[1] + step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - 2(_fd.order - 1),
    ),
    range(
        mx[1] + 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - (_fd.order - 1),
    ),
]
get_mgrid(::tauyz, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] + 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - (_fd.order - 1),
    ),
    range(
        my[1] + 0.5 * step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - (_fd.order - 1),
    ),
    range(
        mx[1] + step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - 2(_fd.order - 1),
    ),
]

get_mgrid(::tauxx, ::FdtdElastic, mz, mx) = [mz, mx]
get_mgrid(::tauzz, ::FdtdElastic, mz, mx) = [mz, mx]
get_mgrid(::tauxz, ::FdtdElastic, mz, mx) = [
    range(
        mz[1] + 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - (_fd.order - 1),
    ),
    range(
        mx[1] + 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - (_fd.order - 1),
    ),
]

get_mgrid(::vx, ::FdtdElastic, mz, my, mx) = [
    mz,
    my,
    range(
        mx[1] - 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) + (_fd.order - 1),
    ),
]
get_mgrid(::vy, ::FdtdElastic, mz, my, mx) = [
    mz,
    range(
        my[1] - 0.5 * step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) + (_fd.order - 1),
    ),
    mx,
]
get_mgrid(::vz, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] - 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) + (_fd.order - 1),
    ),
    my,
    mx,
]

get_mgrid(::vx, ::FdtdElastic, mz, mx) = [
    mz,
    range(
        mx[1] - 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) + (_fd.order - 1),
    ),
]
get_mgrid(::vz, ::FdtdElastic, mz, mx) = [
    range(
        mz[1] - 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) + (_fd.order - 1),
    ),
    mx,
]

get_mgrid(::dvxdx, ::FdtdElastic, mz, my, mx) = [mz, my, mx]
get_mgrid(::dvydy, ::FdtdElastic, mz, my, mx) = [mz, my, mx]
get_mgrid(::dvzdz, ::FdtdElastic, mz, my, mx) = [mz, my, mx]
get_mgrid(::dvxdy, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] + step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - 2(_fd.order - 1),
    ),
    range(
        my[1] + 0.5 * step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - (_fd.order - 1),
    ),
    range(
        mx[1] + 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - (_fd.order - 1),
    ),
]
get_mgrid(::dvxdz, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] + 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - (_fd.order - 1),
    ),
    range(
        my[1] + step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - 2(_fd.order - 1),
    ),
    range(
        mx[1] + 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - (_fd.order - 1),
    ),
]
get_mgrid(::dvydx, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] + step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - 2(_fd.order - 1),
    ),
    range(
        my[1] + 0.5 * step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - (_fd.order - 1),
    ),
    range(
        mx[1] + 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - (_fd.order - 1),
    ),
]
get_mgrid(::dvydz, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] + 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - (_fd.order - 1),
    ),
    range(
        my[1] + 0.5 * step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - (_fd.order - 1),
    ),
    range(
        mx[1] + step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - 2(_fd.order - 1),
    ),
]
get_mgrid(::dvzdx, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] + 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - (_fd.order - 1),
    ),
    range(
        my[1] + step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - 2(_fd.order - 1),
    ),
    range(
        mx[1] + 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - (_fd.order - 1),
    ),
]
get_mgrid(::dvzdy, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] + 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - (_fd.order - 1),
    ),
    range(
        my[1] + 0.5 * step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - (_fd.order - 1),
    ),
    range(
        mx[1] + step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - 2(_fd.order - 1),
    ),
]

get_mgrid(::dvxdx, ::FdtdElastic, mz, mx) = [mz, mx]
get_mgrid(::dvzdz, ::FdtdElastic, mz, mx) = [mz, mx]
get_mgrid(::dvxdz, ::FdtdElastic, mz, mx) = [
    range(
        mz[1] + 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - (_fd.order - 1),
    ),
    range(
        mx[1] + 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - (_fd.order - 1),
    ),
]
get_mgrid(::dvzdx, ::FdtdElastic, mz, mx) = [
    range(
        mz[1] + 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - (_fd.order - 1),
    ),
    range(
        mx[1] + 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - (_fd.order - 1),
    ),
]

get_mgrid(::dtauxxdx, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] + step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - 2(_fd.order - 1),
    ),
    range(
        my[1] + step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - 2(_fd.order - 1),
    ),
    range(
        mx[1] + 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - (_fd.order - 1),
    ),
]
get_mgrid(::dtauyydy, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] + step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - 2(_fd.order - 1),
    ),
    range(
        my[1] + 0.5 * step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - (_fd.order - 1),
    ),
    range(
        mx[1] + step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - 2(_fd.order - 1),
    ),
]
get_mgrid(::dtauzzdz, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] + 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - (_fd.order - 1),
    ),
    range(
        my[1] + step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - 2(_fd.order - 1),
    ),
    range(
        mx[1] + step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - 2(_fd.order - 1),
    ),
]
get_mgrid(::dtauxydx, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] + step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - 2(_fd.order - 1),
    ),
    range(
        my[1] + 0.5 * step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - (_fd.order - 1),
    ),
    range(
        mx[1] + step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - 2(_fd.order - 1),
    ),
]
get_mgrid(::dtauxydy, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] + step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - 2(_fd.order - 1),
    ),
    range(
        my[1] + step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - 2(_fd.order - 1),
    ),
    range(
        mx[1] + 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - (_fd.order - 1),
    ),
]
get_mgrid(::dtauxzdx, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] + 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - (_fd.order - 1),
    ),
    range(
        my[1] + step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - 2(_fd.order - 1),
    ),
    range(
        mx[1] + step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - 2(_fd.order - 1),
    ),
]
get_mgrid(::dtauxzdz, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] + step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - 2(_fd.order - 1),
    ),
    range(
        my[1] + step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - 2(_fd.order - 1),
    ),
    range(
        mx[1] + 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - (_fd.order - 1),
    ),
]
get_mgrid(::dtauyzdy, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] + 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - (_fd.order - 1),
    ),
    range(
        my[1] + step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - 2(_fd.order - 1),
    ),
    range(
        mx[1] + step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - 2(_fd.order - 1),
    ),
]
get_mgrid(::dtauyzdz, ::FdtdElastic, mz, my, mx) = [
    range(
        mz[1] + step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - 2(_fd.order - 1),
    ),
    range(
        my[1] + 0.5 * step(my) * (_fd.order - 1),
        step = step(my),
        length = length(my) - (_fd.order - 1),
    ),
    range(
        mx[1] + step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - 2(_fd.order - 1),
    ),
]

get_mgrid(::dtauxxdx, ::FdtdElastic, mz, mx) = [
    range(
        mz[1] + step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - 2(_fd.order - 1),
    ),
    range(
        mx[1] + 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - (_fd.order - 1),
    ),
]
get_mgrid(::dtauzzdz, ::FdtdElastic, mz, mx) = [
    range(
        mz[1] + 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - (_fd.order - 1),
    ),
    range(
        mx[1] + step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - 2(_fd.order - 1),
    ),
]
get_mgrid(::dtauxzdx, ::FdtdElastic, mz, mx) = [
    range(
        mz[1] + 0.5 * step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - (_fd.order - 1),
    ),
    range(
        mx[1] + step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - 2(_fd.order - 1),
    ),
]
get_mgrid(::dtauxzdz, ::FdtdElastic, mz, mx) = [
    range(
        mz[1] + step(mz) * (_fd.order - 1),
        step = step(mz),
        length = length(mz) - 2(_fd.order - 1),
    ),
    range(
        mx[1] + 0.5 * step(mx) * (_fd.order - 1),
        step = step(mx),
        length = length(mx) - (_fd.order - 1),
    ),
]



