
using GeoPhyInv


# # Intro
#md # ```@docs
#md # AGeom
#md # ```

# # Examples

# Lets create a 2-D `mgrid` for the experiment.
mgrid2 = [range(0.0, stop = 10.0, step = 0.1), range(0.0, stop = 30.0, step = 0.2)];

# Then initialize an acquisition on `mgrid`, where the positions will be randomly chosen.
ageom2 = AGeom(mgrid2, SSrcs(10), Srcs(10), Recs(10));

# Similarly, we can do it for 3D too.
mgrid3 = fill(range(-10, stop=10, step=0.01),3);
ageom3 = AGeom(mgrid3, SSrcs(10), Srcs(10), Recs(10));

# You can check the number of dimensions.
ndims(ageom2), ndims(ageom3) == 2, 3


# For 2D, we can also initialize with one of the predefined acquisitions.
ageom2 = AGeom(mgrid2, :xwell, SSrcs(10), Recs(10));

# By default, the sources and receivers are randomly placed
# on the grid.
# If unsure, test it.
@test (ageom2 ∈ mgrid2)
@test (ageom3 ∈ mgrid3)

# The source and receiver positions can be updated as desired.  
update!(ageom2[1], Srcs(), [0, 1], [10, 20]);
update!(ageom2[1], Recs(), [0, 0], [10, 20]);
update!(ageom2, SSrcs(), [0, 1], [10, 20]);

# It is easy to combine supersources. Now `ageom2` has 20 supersources.
ageom2_new = vcat(ageom2, ageom2);

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

