
function SeisForwExpt(attrib::Symbol; npw=1)
    @assert attrib in [:acou_homo2D, :elastic_homo2D]

    medium = Medium(:elastic_homo2D, 15)
    # update!(medium, [:vp, :rho, :vs], randn_perc = 5)
    # ageom = AGeom(medium.mgrid, :xwell, SSrcs(1), Recs(100))
    ageom = AGeom(medium.mgrid, :xwell, SSrcs(1), Recs(1))

    wav, tgrid = ricker(medium, 10, 1.0)
    rmul!(wav, 1e6)
    srcwav = SrcWav(tgrid, ageom, [:vz])
    update!(srcwav, [:vz], wav)

    tsnaps = tgrid[1:div(length(tgrid), 10):end]
    if (attrib == :acou_homo2D)
        attrib_mod = FdtdAcoustic{FullWave}(:forward, npw) 
    elseif (attrib == :elastic_homo2D)
        attrib_mod = FdtdElastic{FullWave}(:forward, npw)
    end

    return SeisForwExpt(
        attrib_mod,
        tgrid = tgrid,
        tsnaps = tsnaps,
        snaps_field = :vz,
        rfields = [:vz],
        pml_faces = [:xmin, :xmax, :zmin, :zmax],
        rigid_faces = [:xmin, :xmax, :zmin, :zmax],
        verbose = true,
        backprop_flag = :save,
        illum_flag = false,
        ageom = ageom,
        srcwav = srcwav,
        medium = medium,
    )


end

