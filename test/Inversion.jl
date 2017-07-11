# add processors
addprocs(2)
using SIT
using Base.Test

# test function
f4(x::Vector) = (100.0 - x[1])^2 + (50.0 - x[2])^2

x=[100.0, 50.0]
gx1 = similar(x)
gx2 = similar(x)

# computing gradient in parallel
SIT.Inversion.finite_difference!(x->(100.0 - x[1])^2 + (50.0 - x[2])^2, x, gx1, :central)


# computing gradient 
SIT.Inversion.finite_difference!(x->(100.0 - x[1])^2 + (50.0 - x[2])^2, x, gx2, :forward)

@test norm(gx1-[0., 0.]) < 10e-4
@test norm(gx2-[0., 0.]) < 10e-4
@test norm(gx2-gx1) < 10e-4







