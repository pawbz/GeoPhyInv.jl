# This page was generated on DATEOFTODAY
using GeoPhyInv
using Random
using LinearAlgebra
#jl @init_parallel_stencil(3, false, Float32, 4)

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

# Lets initialize records for `:p` field.
srcwav=SrcWav(tgrid, ageom, [:p]);

# Fill the `:P` field of 3rd supersource with random numbers.
randn!(srcwav[3][:p]);

# Often we want to populate the same source wavelet to all
# the supersources and sources.
x=randn(length(tgrid));
update!(srcwav, [:p,], x);

# Populate two different wavelets for first and second supersources.
x1=randn(length(tgrid));
x2=randn(length(tgrid));
update!(srcwav[1], [:p,], x1);
update!(srcwav[2], [:p,], x2);


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

