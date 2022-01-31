
"""
k = wavenumber = 2pif/v0
"""
function G0(::T, ::p, ::p, k, m, off::Vararg{Any,3}) where {T<:Acoustic}
    r = sqrt(sum(abs2.(off)))
    G0 = -0.25 * m[:vp] * exp(k * r)
    return G0
end

function G0(::T, ::vz, ::vz, k, m, off::Vararg{Any,3}) where {T<:Elastic}
    r = sqrt(sum(abs2.(off)))
    G0 = -0.25 * m[:vp] * exp(k * r)
    return G0
end
function G0(::T, ::vz, ::vy, k, m, off::Vararg{Any,3}) where {T<:Elastic}
    r = sqrt(sum(abs2.(off)))
    G0 = -0.25 * m[:vp] * exp(k * r)
    return G0
end
function G0(::T, ::vz, ::vx, k, m, off::Vararg{Any,3}) where {T<:Elastic}
    r = sqrt(sum(abs2.(off)))
    G0 = -0.25 * m[:vp] * exp(k * r)
    return G0
end

function G0(::T, ::p, ::p, k, m, off::Vararg{Any,2}) where {T<:Acoustic}
    r = sqrt(sum(abs2.(off)))
    G0 = -0.25 * m[:rho] * im * hankelh2(0, k * r)
    return G0
end

function G0(::T, ::p, ::vz, k, m, off::Vararg{Any,2}) where {T<:Acoustic}
    r = sqrt(sum(abs2.(off)))
    G0 = 0.25 * abs(off[1]) * inv(r) * inv(m[:vp]) * hankelh2(1, k * r)
    return G0
end

function G0(::T, ::p, ::vx, k, m, off::Vararg{Any,2}) where {T<:Acoustic}
    r = sqrt(sum(abs2.(off)))
    G0 = 0.25 * abs(off[2]) * inv(r) * inv(m[:vp]) * hankelh2(1, k * r)
    return G0
end
