using SIT

# Rosenbrook
f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

x=[0.0,0.0]
gx = similar(x)
SIT.Inversion.finite_difference!(f, x, gx, :central)

println(gx)




