# This page was generated on DATEOFTODAY

using BenchmarkTools
using GeoPhyInv
using Test
using Luxor

# ## Staggered Grid
# Luxor graphics are available to visualize 2-D grids.
# 3-D visualization is not available.
# ### Acoustic
@draw GeoPhyInv.luxor_mgrid(FdtdAcoustic())
# We can avoid some clutter.
@draw GeoPhyInv.luxor_mgrid(FdtdAcoustic(), [:p, :vx, :vz])


# ### Elastic
# Now lets visualize the grids of elastic simulation.
@draw GeoPhyInv.luxor_mgrid(FdtdElastic())
# Grids for only stress and velocity fields.
@draw GeoPhyInv.luxor_mgrid(FdtdElastic(), [:tauxx, :tauxz, :tauzz, :vx, :vz])



# ## Methods

#md # ```@docs
#md # SeisForwExpt 
#md # Base.getindex(::GeoPhyInv.PFdtd, ::Symbol, ::Int)
#md # ```


#md # ```@docs
#md # update!(::GeoPhyInv.PFdtd)
#md # update!(::GeoPhyInv.PFdtd, ::Medium)
#md # update!(::GeoPhyInv.PFdtd, ::SrcWav, ::Any)
#md # ```
#
#
#
#
#
#
#


