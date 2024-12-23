
"""
Modify the input acquisition ageometry 
such that the adjoint source time functions can 
be propagated from the receiver positions.
The number of supersources will remain the same.
All the recievers will be fired as simultaneous sources.
THIS IS NOT USED ANYMORE, delete?
"""
function get_adjoint_ageom(ageom::AGeom)
    return broadcast(ageom) do a
        AGeomss(a.r, a.r, a.nr, a.nr)
    end
end

for dimnames in [zip([:1, :2, :3], [:z, :y, :x]), zip([:1, :2], [:z, :x])]
    grids = broadcast(x -> Symbol(string("m", x)), getindex.(collect(dimnames), 2))
    # get indices for any fields, routed via mgrid adjustment
    @eval function get_proj_matrix(f::Any, attrib_mod, $(grids...), Ps)
        return get_proj_matrix(Ps, get_mgrid(f, attrib_mod, $(grids...)))
    end
end


"""
This function updates spray and interpolation matrices. Call it whenever there is a change in the source positions or 
their fields.
"""
function update!(pass::P_x_worker_x_pw_x_ss, ipw, iss, ageomss::AGeomss, pac)
    update!(pass, ipw, iss, ageomss, pac, Srcs(0))
    update!(pass, ipw, iss, ageomss, pac, Recs(0))
end
function update!(pass::P_x_worker_x_pw_x_ss, ipw, iss, ageomss::AGeomss, pac, ::Srcs)
    @assert ageomss.ns == pac.ageom[ipw][iss].ns

    ssprayw = pass.ssprayw

    sfields = names(pac.srcwav[ipw][iss].d)[1]
    setnames!(ssprayw, sfields, 1)

    for sfield in sfields
        ssprayw[sfield] = get_proj_matrix(eval(sfield)(), pac.attrib_mod, pac.exmedium.grid...,[[s[is] for s in ageomss.s] for is in 1:ageomss.ns])
    end

end
function update!(pass::P_x_worker_x_pw_x_ss, ipw, iss, ageomss::AGeomss, pac, ::Recs)
    @assert ageomss.nr == pac.ageom[ipw][iss].nr

    rinterpolatew = pass.rinterpolatew

    for rfield in pac.rfields
        rinterpolatew[rfield] = get_proj_matrix(eval(rfield)(), pac.attrib_mod, pac.exmedium.grid...,[[r[ir] for r in ageomss.r] for ir in 1:ageomss.nr])
    end
end


# if just one propagating field
update!(pa::PFdtd, ageom::AGeom) = update!(pa, [ageom])
update!(pa::PFdtd, ageom::AGeom, ::Srcs) = update!(pa, [ageom], Srcs(0))
update!(pa::PFdtd, ageom::AGeom, ::Recs) = update!(pa, [ageom], Recs(0))

function update!(pa::PFdtd, ageom::Vector{AGeom})
    update!(pa, ageom, Srcs(0))
    update!(pa, ageom, Recs(0))
end
function update!(pa::PFdtd, ageom::Vector{AGeom}, sr::Srcs)
    for ipw = 1:pa.c.ic[:npw]
        copyto!(pa.c.ageom[ipw], ageom[ipw])
        @sync begin
            for (ip, p) in enumerate(procs(pa.p))
                @async remotecall_wait(p) do
                    pap = localpart(pa.p)[ipw]
                    for is = 1:length(pap.ss)
                        iss = pap.ss[is].iss
                        update!(localpart(pap).ss[is], ipw, iss, ageom[ipw][iss], pa.c, sr)
                    end
                end
            end
        end
    end
end

function update!(pa::PFdtd, ageom::Vector{AGeom}, sr::Recs)
    for ipw = 1:pa.c.ic[:npw]
        copyto!(pa.c.ageom[ipw], ageom[ipw])
        @sync begin
            for (ip, p) in enumerate(procs(pa.p))
                @async remotecall_wait(p) do
                    pap = localpart(pa.p)[ipw]
                    for is = 1:length(pap.ss)
                        iss = pap.ss[is].iss
                        update!(localpart(pap).ss[is], ipw, iss, ageom[ipw][iss], pa.c, sr)
                    end
                end
            end
        end
    end
end

