abstract type Elastic end
abstract type Acoustic end

"""
2-D/ 3-D Linearized forward modeling using a finite-difference simulation of the wave-equation.
"""
struct Born end

"""
Full wavefield modelling; chosen by default.
"""
struct FWmod end

"""
2-D/ 3-D elastic forward modeling using staggered-grid velocity-stress finite-difference formulation.
"""
struct FdtdElastic{T} <: Elastic end
FdtdElastic()=FdtdElastic{FWmod}()
"""
2-D/ 3-D acoustic forward modeling using staggered-grid velocity-stress finite-difference formulation.
"""
struct FdtdAcoustic{T} <: Acoustic end
FdtdAcoustic()=FdtdAcoustic{FWmod}()
"""
2-D/ 3-D viscoacoustic forward modeling using staggered-grid velocity-stress finite-difference formulation.
"""
struct FdtdAcousticVisco end
"""
2-D/ 3-D Linearized forward modeling using a analytical solutions for homogeneous acoustic media.
"""
struct AcousticBorn <: Acoustic end
"""
2-D/ 3-D Linearized forward modeling using a analytical solutions for homogeneous elastic media.
"""
struct ElasticBorn <: Elastic end

# define supersource, source and receiver structs
struct Srcs
    n::Int
end
Srcs() = Srcs(0)
struct SSrcs
    n::Int
end
SSrcs() = SSrcs(0)
struct Recs
    n::Int
end
Recs() = Recs(0)


