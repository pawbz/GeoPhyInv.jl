
function SeisForwExpt(attrib::Symbol)
    @assert attrib in [:acou_homo2D, :elastic_homo2D]

    medium = Medium(:elastic_homo2D, 5)
    update!(medium, [:vp, :rho, :vs], randn_perc = 5)
    ageom = AGeom(medium.mgrid, :xwell, SSrcs(2), Recs(100))

    wav, tgrid = ricker(medium, 10, 1.0)
    srcwav = SrcWav(tgrid, ageom, [:vz])
    update!(srcwav, [:vz], wav)

    tsnaps = tgrid[1:div(length(tgrid), 20):end]
    if (attrib == :acou_homo2D)
        attrib_mod = FdtdAcoustic()
    elseif (attrib == :elastic_homo2D)
        attrib_mod = FdtdElastic()
    end

    return SeisForwExpt(
        attrib_mod,
        tgrid = tgrid,
        tsnaps = tsnaps,
        snaps_field = :vz,
        rfields = [:vx, :vz],
        pml_faces = [:xmin, :xmax, :zmin, :zmax],
        rigid_faces = [:xmin, :xmax, :zmin, :zmax],
        verbose = true,
        backprop_flag = :save,
        illum_flag = false,
        ageom = [ageom],
        srcwav = [srcwav],
        medium = medium,
    )


end

