
using GeoPhyInv
using Test
using SparseArrays


# # Intro
#md # ```@docs
#md # AGeom
#md # ```

# # Examples

# Lets create a `mgrid` for the experiment. 
mgrid = [range(0.0, stop=10.,step=0.1), range(0.0, stop=30.,step=0.2)];

# Then initialize an acquisition on `mgrid`.

ageom=AGeom(mgrid, SSrcs(10), Srcs(10), Recs(10));

# Otherwise, initialize with one of the predefined acquisitions.

ageom=AGeom(mgrid, :xwell, SSrcs(10), Recs(10));

# By default, the sources and receivers are randomly placed
# on the grid.
# If unsure, test it.
@test (ageom ∈ mgrid)

# The source and receiver positions can be updated as desired.  

update!(ageom[1], Srcs(), [0,1], [10,20],);
update!(ageom[1], Recs(), [0,0], [10,20],);
update!(ageom, SSrcs(), [0,1], [10,20], );

# It is easy to combine supersources. Now `ageom2` has 20 supersources.
ageom2=vcat(ageom, ageom);

# # Methods 
#md # ```@docs
#md # update!(::AGeomss, ::Srcs)
#md # update!(::AGeomss, ::Recs)
#md # update!(::AGeom, ::SSrcs)
#md # update!(::AGeom, ::Recs)
#md # Base.in(::AGeom, ::AbstractVector{StepRangeLen})
#md # Base.isequal(::AGeom, ::AGeom)
#md # SparseArrays.SparseMatrixCSC(::AGeomss,::AbstractVector{StepRangeLen})
#md # ```

