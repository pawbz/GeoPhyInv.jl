
"""
```julia
AGeom=Vector{AGeomss,1}
```
An acquisition geometry, bundled into the mutable `AGeom` type,
has to be specified in order to
create any `Expt` variable.
Acquisition has supersources, sources and receivers.
For `AGeom` is a mutable type to define source-receiver geometry. 
Each supersource has `ns` sources that are
injected (or active) simultaneously, and 
a set of `nr` receivers to 
record. 
`AGeomss` is a subtype for each supersource, therefore:
This package provides tools to easily create commonly-used `AGeom` instances, 
after deciding on a spatial grid `mgrid` for the experiment.

## Indexing and Fields
Let `ageom` be an instance of this type, then fields 
can be accessed using:

* `ageom[i]` : acquisition for ith supersource
* `ageom[i].s[p]` : position coordinates of sources where `p ∈ [:z, :x, :y]`
* `ageom[i].r[p]` : position coordinates of receivers
* `ageom[i].ns` : number of sources 
* `ageom[i].nr` : number of receivers


```julia
ageom=AGeom(mgrid, attrib, SSrcs(10), Recs(20))
```
Return pre-defined acquisition geometries (in this case with 10 supersources and 20 receivers), based on `attrib`, on an input mesh `mgrid`.
The sources and receivers positions are same for all supersources, and are not placed on the edges of the mesh.
Choose `attrib::Symbol`
* `=:xwell` cross-well acquisition;
* `=:surf` surface acquisition;
* `=:vsp` vertical seismic profiling;
* `=:rvsp`  reverse vertical seismic profiling;
* `=:downhole` downhole sources and receivers.


TODO: For 3D, just mean(ygrid) is used, need to work on more exotic gallery of 3D acquisitions.
"""
AGeom = Array{AGeomss,1}



function Array{AGeomss,1}(
    mgrid::AbstractArray{T},
    attrib::Symbol,
    ss::SSrcs = SSrcs(1),
    r::Recs = Recs(10),
) where {T<:StepRangeLen}
    N = length(mgrid)
    @assert N ∈ [2, 3]
    # ouput a coordinate given a 1-D grid usign weights a and b
    getp(m, a = 0.5) = a * m[1] + (1 - a) * m[end]
    function getp(m::AbstractArray{T}, a = fill(0.5, length(m))) where {T<:StepRangeLen}
        @assert length(m) == length(a)
        @assert all(a .> 0)
        @assert all(a .< 1)
        return [getp(mgrid[i], a[i]) for i = 1:length(m)]
    end
    ageom = AGeom(mgrid, ss, Srcs(1), r) # just use one source per supersource
    if (attrib == :xwell)
        sp1 = [0.9, 0.9]
        sp2 = [0.1, 0.9]
        rp1 = [0.9, 0.1]
        rp2 = [0.1, 0.1]
    elseif (attrib == :surf)
        sp1 = [0.95, 0.9]
        sp2 = [0.95, 0.1]
        rp1 = [0.95, 0.9]
        rp2 = [0.95, 0.1]
    elseif (attrib == :vsp)
        sp1 = [0.95, 0.9]
        sp2 = [0.95, 0.1]
        rp1 = [0.75, 0.9]
        rp2 = [0.1, 0.9]
    elseif (attrib == :rvsp)
        sp1 = [0.75, 0.9]
        sp2 = [0.1, 0.9]
        rp1 = [0.95, 0.9]
        rp2 = [0.95, 0.1]
    elseif (attrib == :downhole)
        sp1 = [0.75, 0.75]
        sp2 = [0.9, 0.75]
        rp1 = [0.75, 0.75]
        rp2 = [0.9, 0.75]
    elseif (attrib == :microseismic)
        sp1 = [rand(Uniform(0.1,0.2)), rand(Uniform(0.7,0.9))]
        sp2 = [rand(Uniform(0.1,0.2)), rand(Uniform(0.7,0.9))]
        rp1 = [0.95, 0.9]
        rp2 = [0.95, 0.1]
    else
        error("invalid attrib")
    end
    if (N == 2)
        update!(ageom, SSrcs(), getp(mgrid, sp1), getp(mgrid, sp2))
        update!(ageom, Recs(), getp(mgrid, rp1), getp(mgrid, rp2))
    else # just use mean y for 3D 
        update!(
            ageom,
            SSrcs(),
            getp(mgrid, [sp1[1], 0.5, sp1[2]]),
            getp(mgrid, [sp2[1], 0.5, sp2[2]]),
        )
        update!(
            ageom,
            Recs(),
            getp(mgrid, [rp1[1], 0.5, rp1[2]]),
            getp(mgrid, [rp2[1], 0.5, rp2[2]]),
        )
    end
    return ageom
end

