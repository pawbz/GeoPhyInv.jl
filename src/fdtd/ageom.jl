
"""
Modify the input acquisition ageometry 
such that the adjoint source time functions can 
be propagated from the receiver positions.
The number of supersources will remain the same.
All the recievers will be fired as simultaneous sources.
"""
function get_adjoint_ageom(ageom::AGeom)
    return broadcast(ageom) do a
        AGeomss(a.r, a.r, a.nr, a.nr)
    end
end



"""
This function updates spray and interpolation matrices. Call it whenever there is a change in the source positions or 
their fields.
"""
function update!(pass::P_x_worker_x_pw_x_ss, ipw, iss, ageomss::AGeomss, pac)
    update!(pass, ipw, iss, ageomss, pac, Srcs())
    update!(pass, ipw, iss, ageomss, pac, Recs())
end
function update!(pass::P_x_worker_x_pw_x_ss, ipw, iss, ageomss::AGeomss, pac, ::Srcs)
    @assert ageomss.ns == pac.ageom[ipw][iss].ns

    ssprayw = pass.ssprayw

    sfields = names(pac.srcwav[ipw][iss].d)[1]
    setnames!(ssprayw, sfields, 1)

    for sfield in sfields
        I = Vector{Int64}()
        J = Vector{Int64}()
        V = Vector{Data.Number}()
        for is = 1:ageomss.ns
            In, Jn, Ln, Vn = findnz(
                pac,
                Srcs(),
                eval(sfield)(),
                pac.attrib_mod,
                pac.exmedium.mgrid...,
                [s[is] for s in ageomss.s],
            )
            I = vcat(I, Ln)
            J = vcat(J, fill(is, length(Ln)))
            V = vcat(V, Data.Number.(Vn))
        end
        M = prod(length.(get_mgrid(eval(sfield)(), pac.attrib_mod, pac.exmedium.mgrid...)))
        if (_fd.use_gpu)
            ssprayw[sfield] = CuSparseMatrixCSC(sparse(I, J, V, M, ageomss.ns))
        else
            ssprayw[sfield] = sparse(I, J, V, M, ageomss.ns)
        end
    end

end
function update!(pass::P_x_worker_x_pw_x_ss, ipw, iss, ageomss::AGeomss, pac, ::Recs)
    @assert ageomss.nr == pac.ageom[ipw][iss].nr

    rinterpolatew = pass.rinterpolatew

    for rfield in pac.rfields
        I = Vector{Int64}()
        J = Vector{Int64}()
        V = Vector{Data.Number}()
        for ir = 1:ageomss.nr
            In, Jn, Ln, Vn = findnz(
                pac,
                Recs(),
                eval(rfield)(),
                pac.attrib_mod,
                pac.exmedium.mgrid...,
                [r[ir] for r in ageomss.r],
            )
            J = vcat(J, Ln)
            I = vcat(I, fill(ir, length(Ln)))
            V = vcat(V, Data.Number.(Vn))
        end
        M = prod(length.(get_mgrid(eval(rfield)(), pac.attrib_mod, pac.exmedium.mgrid...)))
        if (_fd.use_gpu)
            rinterpolatew[rfield] = CuSparseMatrixCSC(sparse(I, J, V, ageomss.nr, M))
        else
            rinterpolatew[rfield] = sparse(I, J, V, ageomss.nr, M)
        end
    end

end


# if just one propagating field
update!(pa::PFdtd, ageom::AGeom) = update!(pa, [ageom])
update!(pa::PFdtd, ageom::AGeom, ::Srcs) = update!(pa, [ageom], Srcs())
update!(pa::PFdtd, ageom::AGeom, ::Recs) = update!(pa, [ageom], Recs())

function update!(pa::PFdtd, ageom::Vector{AGeom})
    update!(pa, ageom, Srcs())
    update!(pa, ageom, Recs())
end
function update!(pa::PFdtd, ageom::Vector{AGeom}, ::Srcs)
    for ipw = 1:pa.c.ic[:npw]
        copyto!(pa.c.ageom[ipw], ageom[ipw])
        @sync begin
            for (ip, p) in enumerate(procs(pa.p))
                @async remotecall_wait(p) do
                    pap = localpart(pa.p)[ipw]
                    for is = 1:length(pap.ss)
                        iss = pap.ss[is].iss
                        update!(localpart(pap).ss[is], ipw, iss, ageom[ipw][iss], pa.c)
                    end
                end
            end
        end
    end
end

function update!(pa::PFdtd, ageom::Vector{AGeom}, ::Recs)
    for ipw = 1:pa.c.ic[:npw]
        copyto!(pa.c.ageom[ipw], ageom[ipw])
        @sync begin
            for (ip, p) in enumerate(procs(pa.p))
                @async remotecall_wait(p) do
                    pap = localpart(pa.p)[ipw]
                    for is = 1:length(pap.ss)
                        iss = pap.ss[is].iss
                        update!(localpart(pap).ss[is], ipw, iss, ageom[ipw][iss], pa.c)
                    end
                end
            end
        end
    end
end

