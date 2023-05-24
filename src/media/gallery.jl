function update_using_δ(medium, δ)
    if (δ == 0.0)
        return medium
    elseif (δ > 0.0)
        mgrid_out = broadcast(x -> range(x[1], stop=x[end], step=δ), medium.grid)
        medium_out = update(medium, mgrid_out)
        return medium_out
    else
        error("invalid δ")
    end
end


vp0(T, mgrid) = vp(fill(T(2500.0), length.(mgrid)...))
rho0(T, mgrid) = rho(fill(T(2500.0), length.(mgrid)...))
vs0(T, mgrid) = vs(fill(T(1500.0), length.(mgrid)...))
function AcousticMedium{T,N}(::Homogeneous, δ=0.0) where {T,N}
    mgrid = repeat([range(-1000.0, stop=1000.0, length=201)], N)
    medium = AcousticMedium(mgrid, vp0(T, mgrid), rho0(T, mgrid))
    return update_using_δ(medium, δ)
end
function ElasticMedium{T,N}(::Homogeneous, δ=0.0) where {T,N}
    mgrid = repeat([range(-1000.0, stop=1000.0, length=201)], N)
    medium = ElasticMedium(mgrid, vp0(T, mgrid), vs0(T, mgrid), rho0(T, mgrid))
    return update_using_δ(medium, δ)
end
# defaults to Float32, and 2D
AcousticMedium(a::Gallery, δ=0.0) = AcousticMedium{Data.Number,2}(a, δ)
ElasticMedium(a::Gallery, δ=0.0) = ElasticMedium{Data.Number,2}(a, δ)


# elseif ((attrib == :pizza))
#     mgrid = repeat([range(-1000.0, stop = 1000.0, length = 201)], 2)
#     model = Medium(mgrid, [:vp, :rho])
#     update!(model, [:vp, :rho], [vp0, rho0])
#     fill!(model)
#     # add perturbations
#     for rect_loc in [[500.0, 0.0], [0.0, 500.0], [-500.0, 0.0], [0.0, -500.0]]
#         update!(
#             model,
#             [:vp, :rho],
#             rectangle = [rect_loc, rect_loc .+ 100.0],
#             perc = 20.0,
#         )
#     end

function get_marmousi()
    c = h5open(joinpath(marmousi_folder, "marmousi2.h5"), "r") do file
        vpm = read(file, "vp")
        vsm = read(file, "vs")
        rhom = read(file, "rho")
        xgrid = read(file, "xgrid")
        zgrid = read(file, "zgrid")
        mgrid = [
            range(zgrid[1], stop=zgrid[end], length=size(vpm, 1)),
            range(xgrid[1], stop=xgrid[end], length=size(vpm, 2)),
        ]
    return mgrid, vpm, vsm, rhom
    end
end

function AcousticMedium{T,N}(::Marmousi2, δ=0.0) where {T,N}
    mgrid, vpm, vsm, rhom = get_marmousi()
    return AcousticMedium(mgrid, vp(T.(vpm)), rho(T.(rhom .* 1000)))
end
function ElasticMedium{T,N}(::Marmousi2, δ=0.0) where {T,N}
    mgrid, vpm, vsm, rhom = get_marmousi()
    return ElasticMedium(mgrid, vp(T.(vpm)), vs(T.(vsm)), rho(T.(rhom .* 1000)))
end


# elseif (attrib == :marmousi2_small)
#     fmodel = Medium(:marmousi2)
#     fmgrid = fmodel.grid
#     mgrid = [
#         range(500, stop = 3500, step = step(fmgrid[1])),
#         range(4000, step = step(fmgrid[2]), stop = 13000),
#     ]
#     model = update(fmodel, mgrid)

# elseif (attrib == :overthrust)
#     vp = []
#     for i in range(1, 8, step = 1)
#         file_ =
#             h5open(joinpath(overthrust_folder, "overthrust_" * string(i) * ".h5"), "r")
#         vp_sect = read(file_, "vp")
#         if i == 1
#             vp = deepcopy(vp_sect)
#         else
#             vp = cat(vp, vp_sect, dims = 3)
#         end
#         close(file_)
#     end
#     file_ = h5open(joinpath(overthrust_folder, "overthrust_mgrid.h5"), "r")
#     xgrid = read(file_, "xgrid")
#     ygrid = read(file_, "ygrid")
#     zgrid = read(file_, "zgrid")
#     close(file_)
#     mgrid = [
#         range(zgrid[1], stop = zgrid[end], length = size(vp, 1)),
#         range(ygrid[1], stop = ygrid[end], length = size(vp, 2)),
#         range(xgrid[1], stop = xgrid[end], length = size(vp, 3)),
#     ]
#     model = Medium(mgrid, [:vp, :rho])
#     rho0 = [1500.0, 2500.0]
#     update!(model, [:rho], [rho0])
#     fill!(model)
#     copyto!(model[:vp], vp)
#     update!(model, bfrac)

