
"""
2-D/ 3-D elastic forward modeling using staggered-grid velocity-stress finite-difference formulation.
"""
struct FdtdElastic end
"""
2-D/ 3-D acoustic forward modeling using staggered-grid velocity-stress finite-difference formulation.
"""
struct FdtdAcoustic end
"""
2-D/ 3-D viscoacoustic forward modeling using staggered-grid velocity-stress finite-difference formulation.
"""
struct FdtdAcousticVisco end
"""
2-D Linearized forward modeling using a finite-difference simulation of the acoustic wave-equation.
"""
struct FdtdAcousticBorn end


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


