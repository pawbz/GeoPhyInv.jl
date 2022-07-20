using BenchmarkTools
using GeoPhyInv
using Test
using Luxor

@draw GeoPhyInv.luxor_mgrid(FdtdAcoustic())

@draw GeoPhyInv.luxor_mgrid(FdtdAcoustic(), [:p, :vx, :vz])

@draw GeoPhyInv.luxor_mgrid(FdtdElastic())

@draw GeoPhyInv.luxor_mgrid(FdtdElastic(), [:tauxx, :tauxz, :tauzz, :vx, :vz])

