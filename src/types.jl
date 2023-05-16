abstract type Elastic end
abstract type Acoustic end

"""
2-D/ 3-D Linearized forward modeling using a finite-difference simulation of the wave-equation.
"""
struct Born end

"""
Full wavefield modelling; chosen by default.
"""
struct FullWave end

"""
2-D/ 3-D elastic forward modeling using staggered-grid velocity-stress finite-difference formulation.
"""
mutable struct FdtdElastic{T} <: Elastic
    mode::Symbol
    npw::Int
end
# FdtdElastic()=FdtdElastic{FullWave}()
"""
2-D/ 3-D acoustic forward modeling using staggered-grid velocity-stress finite-difference formulation.
"""
mutable struct FdtdAcoustic{T} <: Acoustic
    mode::Symbol
    npw::Int
end

# by default, only one propagating wavefield is choosen
# mode, either :forward, :forward_save or :adjoint, is choosen 
for k in [:FdtdAcoustic, :FdtdElastic]
    @eval function $k()
        return $k{FullWave}(:forward, 1)
    end
    @eval function $k{T}() where {T<:FullWave}
        return $k{T}(:forward, 1)
    end
    @eval function $k{T}() where {T<:Born}
        return $k{T}(:forward, 2)
    end
    @eval function $k(mode)
        if (mode == :forward)
            return $k{FullWave}(:forward, 1)
        else
            return $k{FullWave}(mode, 2)
        end
    end
end

# FdtdAcoustic()=FdtdAcoustic{FullWave}()
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


