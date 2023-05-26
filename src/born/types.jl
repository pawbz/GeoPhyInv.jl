


mutable struct PBorn{T}
    attrib_mod::T
    medium::Medium
    ageom::AGeom
    srcwav::Srcs
    fc::NamedStack{Float64}
    ic::NamedStack{Int64}
    cc::NamedStack{ComplexF64}
    fnpow2grid::FFTW.Frequencies{Float64}
    w::NamedStack{Vector{ComplexF64}}
    m::Vector{NamedStack{ComplexF64}}
    data::Recs
end