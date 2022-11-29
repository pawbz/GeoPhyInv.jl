

# """
# get source spray weights
# """
# function get_spray_weights!(
#     weights::AbstractArray, #denomI::AbstractArray, 
#     ix1::AbstractArray{Int64},
#     ix2::AbstractArray{Int64},
#     iz1::AbstractArray{Int64},
#     iz2::AbstractArray{Int64},
#     mesh_x::AbstractArray,
#     mesh_z::AbstractArray,
#     xval::Float64,
#     zval::Float64,
# )


#     ix1[1], ix2[1] = indminn(mesh_x, xval, 2)
#     iz1[1], iz2[1] = indminn(mesh_z, zval, 2)

#     denomI = inv((mesh_x[ix2[1]] - mesh_x[ix1[1]]) * (mesh_z[iz2[1]] - mesh_z[iz1[1]]))
#     "for iz1, ix1"
#     weights[1] = (mesh_x[ix2[1]] - xval) * (mesh_z[iz2[1]] - zval) * denomI
#     "for iz1, ix2"
#     weights[2] = (xval - mesh_x[ix1[1]]) * (mesh_z[iz2[1]] - zval) * denomI
#     "for iz2, ix1"
#     weights[3] = (mesh_x[ix2[1]] - xval) * (zval - mesh_z[iz1[1]]) * denomI
#     "for iz2, ix2"
#     weights[4] = (xval - mesh_x[ix1[1]]) * (zval - mesh_z[iz1[1]]) * denomI

#     return weights, denomI, ix1, ix2, iz1, iz2
# end

