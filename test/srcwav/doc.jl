
using GeoPhyInv
using LinearAlgebra
using Random

# # Intro

#md # ```@docs
#md # SrcWav
#md # ```


# # Examples
# Define an acquisition geometry.
mgrid = [range(0.0, stop=10.,step=0.1), range(0.0, stop=30.,step=0.2)];
ageom=AGeom(mgrid, SSrcs(10), Srcs(10), Recs(10));


# Need a time grid.
tgrid=range(0, stop=1.0, step=0.1);

# Lets initialize records for `:P` field.
srcwav=SrcWav(tgrid, ageom, [:P]);

# Fill the `:P` field of 3rd supersource with random numbers.
Random.randn!(srcwav[3][:P]);

# Often we want to populate the same source wavelet to all
# the supersources and sources.
x=randn(length(tgrid));
update!(srcwav, [:P,], x);

# Populate two different wavelets for first and second supersources.
x1=randn(length(tgrid));
x2=randn(length(tgrid));
update!(srcwav[1], [:P,], x1);
update!(srcwav[2], [:P,], x2);


# Scale `srcwav` by a scalar overwriting it in-place.
rmul!(srcwav, 2.0);


# # Methods 
# Most of the methods listed below are also applicable to individual elements of `srcwav`.

#md # ```@docs
#md # update!(GeoPhyInv.VNamedD, ::Vector{Symbol}, ::AbstractArray)
#md # Base.reverse!(::GeoPhyInv.VNamedD)
#md # Base.iszero(::GeoPhyInv.VNamedD)
#md # Base.isequal(::GeoPhyInv.VNamedD, ::GeoPhyInv.VNamedD)
#md # GeoPhyInv.issimilar(::GeoPhyInv.VNamedD, ::GeoPhyInv.VNamedD)
#md # Base.vec(::GeoPhyInv.VNamedD)
#md # Random.randn!(::GeoPhyInv.VNamedD)
#md # Base.fill!(::GeoPhyInv.VNamedD, ::Float64)
#md # LinearAlgebra.dot(::GeoPhyInv.VNamedD, ::GeoPhyInv.VNamedD)
#md # LinearAlgebra.rmul!(::GeoPhyInv.VNamedD, ::Number)
#md # Base.copyto!(::GeoPhyInv.VNamedD, ::GeoPhyInv.VNamedD)
#md # Base.copyto!(::AbstractVector{Float64}, ::GeoPhyInv.VNamedD)
#md # Base.copyto!(::GeoPhyInv.VNamedD, ::AbstractVector{Float64})
#md # ```

