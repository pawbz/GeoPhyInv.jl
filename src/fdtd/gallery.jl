
function SeisForwExpt(attrib_mod::FdtdAcoustic, ::Homogeneous; snaps_field=nothing, randn_perc=0.0)

    medium = AcousticMedium(Homogeneous(), 15.0)
    update!(medium, [:vp, :rho], randn_perc=randn_perc)

    ageom = AGeom(medium.grid, :xwell, SSrcs(1), Recs(100))
    # ageom = AGeom(medium.grid, :xwell, SSrcs(1), Recs(1))
    tgrid = range(0, stop=2.0, step=0.003)
    wav = ricker(10.0, tgrid)
    rmul!(wav, 1e6)
    srcwav = Srcs(tgrid, ageom, [:vz])
    update!(srcwav, [:vz], wav)

    tsnaps = tgrid
    return SeisForwExpt(
        attrib_mod,
        tgrid=tgrid,
        tsnaps=tsnaps,
        snaps_field=snaps_field,
        rfields=[:vz],
        pml_faces=[:xmin, :xmax, :zmin, :zmax],
        rigid_faces=[:xmin, :xmax, :zmin, :zmax],
        verbose=true,
        backprop_flag=:save,
        illum_flag=false,
        ageom=ageom,
        srcwav=srcwav,
        medium=medium,
    )
end
SeisForwExpt(attrib_mod::Union{FdtdAcoustic,FdtdElastic}, ::RandScatterer; snaps_field=nothing) = SeisForwExpt(attrib_mod, Homogeneous(); snaps_field=snaps_field, randn_perc=5)
#     @assert attrib in [:acou_homo2D, :elastic_homo2D]


