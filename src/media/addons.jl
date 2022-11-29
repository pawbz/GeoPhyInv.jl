# perturb an array m by percentage of m0, when mask is true
function perturb!(m, m0, perc, mask::Function)
    for i in CartesianIndices(m)
        ii = Tuple(i)
        if (mask(ii))
            m[i] = m[i] + perc * 1e-2 * m0
        end
    end
end
function perturb!(m, m0, perc, idx::CartesianIndices)
    for i in idx
        m[i] = m[i] + perc * 1e-2 * m0
    end
end


"""
In-place method to add a rectangular box or cuboid in the medium. Its position is specified by its corners  
and perc denotes the percentage of the perturbation.

```julia
update!(medium, field; rectangle=[[0,10],[10,10]], perc=10) # for 2D medium
update!(medium, field; rectangle=[[0,10,15],[10,10,15]], perc=10) # for 3D medium
```

In-place method to add random noise to medium.

```julia
update!(medium, fields; randn_perc=10)
```
"""
function update!(
    medium::Medium,
    fields::Vector{Symbol};
    rectangle = nothing,
    perc = 0.0,
    randn_perc = 0.0,
)

    mgrid = medium.mgrid
    if (!(rectangle === nothing))
        @assert length(rectangle)==2 
        @assert length.(rectangle)==fill(ndims(medium),2)
        rectangle_indices = CartesianIndices(
            Tuple([
                argmin(abs.(mg .- rectangle[1][im])):argmin(abs.(mg .- rectangle[2][im])) for (im, mg) in enumerate(mgrid)
            ]),
        )
        if (!iszero(perc))
            for field in fields
                m = medium.m[field]
                m0 = medium.ref[field]
                perturb!(m, m0, perc, rectangle_indices)
            end
        end
    end
    if (!iszero(randn_perc))
        for field in fields
            m = medium.m[field]
            m0 = medium.ref[field]
            for i in eachindex(m)
                m[i] = χ(Float64(χ(m[i], m0, 1) + randn() * randn_perc * 1e-2), m0, -1)
            end

        end
    end

    return medium
end





# # Keyword Arguments

# * `point_loc::Vector{Float64}=[0., 0.,]` : approx location of point pert
# * `point_pert::Float64=0.0` : perturbation at the point scatterer
# * `ellip_loc::Vector{Float64}=nothing` : location of center of perturbation, [z, x]
# * `ellip_rad::Float64=0.0` : size of elliptic perturbation
# * `ellip_pert::Float64=0.1` : perturbation inside the ellipse
# * `ellip_α=0.0` : rotate the ellipse
# * `rect_loc::Array{Float64}=nothing` : rectangle location, [zmin, xmin, zmax, xmax]
# * `rect_pert::Float64=0.1` : perturbation in a rectangle
# * `constant_pert::Float64=0.0` : constant perturbation 
# * `randn_pert::Float64=0.0` : percentage of reference values for additive random noise
# * `fields::Vector{Symbol}=[:χvp,:χrho,:χvs]` : which fields are to be modified?
# * `onlyin` : `medium` is modified only when field values are in these ranges 
# """
# function update!(
#     medium::Medium,
#     fields::Vector{Symbol};
#     point_loc = fill(0.0, ndims(medium)),
#     point_pert::Float64 = 0.0,
#     ellip_loc = fill(0.0, ndims(medium)),
#     ellip_rad = 0.0,
#     ellip_pert::Real = 0.0,
#     ellip_α = 0.0,
#     rect_loc = fill(0.0, ndims(medium)),
#     rect_pert::Float64 = 0.0,
#     constant_pert::Float64 = 0.0,
#     randn_perc::Real = 0.0,
#     layerlocations = nothing,
# )

#     # only editing basic fields is allowed
#     for field in fields
#         @assert field ∈ [:vp, :vs, :rho, :Q]
#     end
#     rect_loc = convert.(Float64, rect_loc)
#     ellip_loc = convert.(Float64, ellip_loc)

#     temp = zeros(length.(medium.mgrid)...)

#     ipointloc = [
#         Interpolation.indminn(medium.mgrid[i], Float64(point_loc[i]), 1)[1] for
#         i = 1:length(medium.mgrid)
#     ]
#     temp[ipointloc...] += point_pert

#     if (!(ellip_pert == 0.0))
#         α = ellip_α * pi / 180.0
#         # circle or ellipse
#         rads =
#             (length(ellip_rad) == 1) ? [ellip_rad[1], ellip_rad[1]] :
#             [ellip_rad[1], ellip_rad[2]]

#         temp += [
#             (
#                 (
#                     (
#                         (
#                             (medium.mgrid[2][ix] - ellip_loc[2]) * cos(α) +
#                             (medium.mgrid[1][iz] - ellip_loc[1]) * sin(α)
#                         )^2 * inv(rads[1]^2) +
#                         (
#                             (-medium.mgrid[1][iz] + ellip_loc[1]) * cos(α) +
#                             (medium.mgrid[2][ix] - ellip_loc[2]) * sin(α)
#                         )^2 * inv(rads[2]^2)
#                     ) <= 1.0
#                 ) ? Float64(ellip_pert) : 0.0
#             ) for iz = 1:length(medium.mgrid[1]), ix = 1:length(medium.mgrid[2])
#         ]
#     end
#     if (!(rect_pert == 0.0))
#         temp += [
#             (
#                 (
#                     (medium.mgrid[2][ix] - rect_loc[4]) * (medium.mgrid[2][ix] - rect_loc[2]) < 0.0
#                 ) & (
#                     (medium.mgrid[1][iz] - rect_loc[3]) * (medium.mgrid[1][iz] - rect_loc[1]) < 0.0
#                 )
#             ) ? rect_pert : 0.0 for iz = 1:length(medium.mgrid[1]),
#             ix = 1:length(medium.mgrid[2])
#         ]
#     end
#     if (!(constant_pert == 0.0))
#         temp .+= constant_pert
#     end


#     for (iff, field) in enumerate(fields)
#         m = medium.m[field]
#         if (!(layerlocations === nothing))
#             onlyatvalues = []
#             for ipos in layerlocations
#                 ipx = Interpolation.indminn(medium.mgrid[2], Float64(ipos[2]), 1)[1]
#                 ipz = Interpolation.indminn(medium.mgrid[1], Float64(ipos[1]), 1)[1]
#                 push!(onlyatvalues, m[ipz, ipx])
#             end
#         end
#         for i in eachindex(m)
#             # if(((m[i]-onlyin[iff][1])*(m[i]-onlyin[iff][2]))<0.0)
#             if (layerlocations === nothing)
#                 m[i] += temp[i]
#             else
#                 if (m[i] ∈ onlyatvalues)
#                     m[i] += temp[i]
#                 end
#             end
#         end
#     end

#     return nothing
# end