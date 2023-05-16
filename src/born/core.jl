
include("types.jl")
include("homo.jl")
include("getprop.jl")


SeisForwExpt(attrib_mod::Union{AcousticBorn,ElasticBorn}, args1...; args2...) =
    PBorn(attrib_mod, args1...; args2...)

function PBorn(
    attrib_mod;
    medium::Medium = nothing,
    Q = nothing,
    medium_pert::Medium = medium,
    scatter_flag::Bool = false,
    tgrid::StepRangeLen = nothing,
    ageom::AGeom = nothing,
    srcwav::SrcWav = nothing,
    rfields::Vector{Symbol} = [:vz],
    sflag::Int64 = 2,
)
    issimilar(ageom, srcwav) || error("ageom and srcwav mismatch")

    nt = (length(tgrid) == length(srcwav[1].grid)) ? length(tgrid) : error("srcwav tgrid")
    np2 = nextpow(2, 2 * nt)
    fnpow2grid = FFTW.fftfreq(np2, inv(step(tgrid)))

    if (:Q ∈ names(medium.m)[1])
        vp0 = medium.ref[:vp]
        rho0 = medium.ref[:rho]
        Ku = rho0 * vp0 * vp0
        tau_epsilon = medium[:tau_epsilon]
        tau_sigma = medium[:tau_sigma]

        Kc = complexK(Ku, abs.(2.0 * pi .* fnpow2grid), tau_sigma, tau_epsilon)

        rho0 = fill(rho0, length(fnpow2grid))
        vp0 = sqrt.(Kc .* inv.(rho0))
    else
        m = [
            NamedArray(complex.([medium.ref[:vp], medium.ref[:rho]]), [:vp, :rho]) for
            i in fnpow2grid
        ]
    end


    cc = NamedArray([complex(0.0, 0.0)], [:dummy])

    if (scatter_flag)
        mesh_x = medium_pert.mgrid[2]
        mesh_z = medium_pert.mgrid[1]
        nz = length(mesh_z)
        nx = length(mesh_x)
        δx = step(mesh_x)
        δz = step(mesh_z)

        δmodKI = medium_pert[:invK] - (vp0 * vp0 * rho0)^(-1)
        δmodrr = medium_pert[:rhoI] - (rho0)^(-1)
    end


    ic = NamedArray(
        vcat(length.(medium.mgrid), [length(tgrid), np2]),
        vcat(dim_names(ndims(medium), "n"), [:nt, :np2]),
    )

    fc = NamedArray([step(tgrid)], [:dt])


    data = Records(tgrid, ageom, rfields)
    w = NamedArray(vcat([complex.(zeros(np2), zeros(np2)) for i = 1:3]), [:G, :S, :D])

    return PBorn(attrib_mod, medium, ageom, srcwav, fc, ic, cc, fnpow2grid, w, m, data)

end

function update!(pa::PBorn)
    ageom = pa.ageom
    nss = length(ageom)
    data = pa.data
    m = pa.m
    D = pa.w[:D]
    S = pa.w[:S]
    G = pa.w[:G]

    for iss = 1:nss
        for rfield in names(data[iss].d)[1]
            d = data[iss].d[rfield]
            for ir = 1:ageom[iss].nr
                fill!(D, complex(0.0))
                for is = 1:ageom[iss].ns
					off=offset(ageom[iss], is, ir)
					(sum(abs2.(off))==0) && error("distance cannot be zero")
                    sfield = names(pa.srcwav[iss].d)[1][1] # only one field for the source?
                    s = pa.srcwav[iss].d[sfield]
                    fill!(S, complex(0.0))
                    fill!(G, complex(0.0))
                    # zero pad wavelet
                    for it = 1:pa.ic[:nt]
                        S[it] = complex(s[it, is])
                    end
                    FFTW.fft!(S) # source wavelet in the frequency domain
                    # analytical expression for every frequency, except zero
                    G[1] = complex.(0.0, 0.0)
                    for iω = 2:pa.ic[:np2]
                        ω = 2.0 * pi * abs(pa.fnpow2grid[iω])
                        k = ω * inv(m[iω][:vp])
    					(k == 0.0) && error("wavenumber cannot be zero")
                        # if (scatter_flag)
                        #     term = complex(0.0, 0.0)
                        #     for ix = 1:nx
                        #         @simd for iz = 1:nz
                        #             if (δmodKI[iz, ix] ≠ 0.0)
                        #                 term +=
                        #                     (
                        #                         G0_homo_acou(
                        #                             sx[is] - mesh_x[ix],
                        #                             sz[is] - mesh_z[iz],
                        #                             k,
                        #                             rho0,
                        #                         )[1] * G0_homo_acou(
                        #                             rx[ir] - mesh_x[ix],
                        #                             rz[ir] - mesh_z[iz],
                        #                             k,
                        #                             rho0,
                        #                         )[1] .* ω .* ω .* δmodKI[iz, ix]
                        #                     ) *
                        #                     δx *
                        #                     δz
                        #                 # factor due to integration
                        #             end
                        #         end
                        #     end

                        # else
                        # if (src_flag == 2)
                        term = G0(
                            pa.attrib_mod,
                            eval(sfield)(),
                            eval(rfield)(),
                            k,
                            m[iω],
                            off...,
                        )
                        # elseif (src_flag == 1)
                        # G[iω] =
                        # G0_homo_acou(x, z, k, rho0[iω]) * im * abs(pa.fnpow2grid[iω])
                        # else
                        # error("invalid src_flag")
                        # end
                        # end
                        if (pa.fnpow2grid[iω] > 0)
                            G[iω] = term
                        elseif (pa.fnpow2grid[iω] < 0)
                            G[iω] = conj(term)
                        end

                    end
                    # convolution and stack over simultaneous sources
                    D += G .* S
                end

                # back to time domain
                FFTW.ifft!(D)

                # truncate
                for it = 1:pa.ic[:nt]
                    d[it, ir] = real(D[it])
                end
            end
        end
    end
    return data
end

