"""
## TODO:
* add dimension checks to interp_spray!
Reference: https://www.ibiblio.org/e-notes/Splines/bezier.html
"""
module Interpolation


using LinearMaps
using LinearAlgebra
using SparseArrays



abstract type Kernel end

mutable struct Kernel_2D{T<:Real} <: Kernel
    F1::SparseMatrixCSC{T,Int64}
    F2::SparseMatrixCSC{T,Int64}
    xx::Matrix{T}
    X::Matrix{T}
    Xi::Matrix{T}
end

mutable struct Kernel_1D{T<:Real} <: Kernel
    F::SparseMatrixCSC{T,Int64}
end


function Kernel(x::T, xi::T, Battrib::Symbol=:B1) where {T}



    # constructor for y=Ax
    function interp!(y, x, pa)
        for i in eachindex(y)
            y[i] = 0.0
        end
        Interpolation.interp_spray!(x, y, pa, :interp)
    end

    if (length(x) == 1)
        pa = P_core(x, xi, Battrib)
        F1 = LinearMap((y, x) -> interp!(y, x, pa),
            nothing, length(xi[1]), length(x[1]), ismutating=true)
        F1 = SparseArrays.sparse(F1)
        return Kernel_1D(F1)
    elseif (length(x) == 2)
        pa1 = P_core(x[2:2], xi[2:2], Battrib)
        F1 = LinearMap((y, x) -> interp!(y, x, pa1),
            nothing, length(xi[2]), length(x[2]), ismutating=true)
        F1 = SparseArrays.sparse(F1)
        pa2 = P_core(x[1:1], xi[1:1], Battrib)
        F2 = LinearMap((y, x) -> interp!(y, x, pa2),
            nothing, length(xi[1]), length(x[1]), ismutating=true)
        F2 = SparseArrays.sparse(F2)
        xx = zeros(size(F1, 1), length(x[1]))
        # temporary storage matrices
        Xi = zeros(size(F1, 1), size(F2, 1))
        X = zeros(length(x[2]), length(x[1]))
        return Kernel_2D(F1, F2, xx, X, Xi)
    else
        error("dimension should be <=2")

    end

end

function interp_spray!(y, yi, pa::Kernel_1D, attrib)
    if (attrib == :interp)
        mul!(yi, pa.F, y)
    elseif (attrib == :spray)
        mul!(y, adjoint(pa.F), yi)
    end
end
function interp_spray!(y::Matrix{T}, yi::Matrix{T}, pa::Kernel_2D{T}, attrib) where {T<:Real}
    if (attrib == :interp)
        mul!(pa.xx, pa.F1, y)
        mul!(yi, pa.xx, adjoint(pa.F2))
    elseif (attrib == :spray)
        mul!(pa.xx, yi, pa.F2)
        mul!(y, adjoint(pa.F1), pa.xx)
    end
end


"""
Method when the 2D arrays that are to interpolated are in the form of vectors.
Note that multiple 2D arrays can be vcat'ed in these vectors, given by nmod.
"""
function interp_spray!(y::Vector{T}, yi::Vector{T}, pa::Kernel_2D{T}, attrib, nmod::Int) where {T<:Real}
    n = div(length(y), nmod)
    ni = div(length(yi), nmod)
    itot = 0
    if (attrib == :interp)
        for i in 1:nmod
            for ii in eachindex(pa.X)
                pa.X[ii] = y[ii+(i-1)*n]
            end
            interp_spray!(pa.X, pa.Xi, pa, attrib)
            for ii in eachindex(pa.Xi)
                yi[ii+(i-1)*ni] = pa.Xi[ii]
            end
        end
    elseif (attrib == :spray)
        for i in 1:nmod
            for ii in eachindex(pa.Xi)
                pa.Xi[ii] = yi[ii+(i-1)*ni]
            end
            interp_spray!(pa.X, pa.Xi, pa, attrib)
            for ii in eachindex(pa.X)
                y[ii+(i-1)*n] = pa.X[ii]
            end
        end

    end
end


include("misc.jl")
include("core.jl")
include("weights.jl")



end # module
