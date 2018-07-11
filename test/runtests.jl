addprocs(4)
using JuMIT
using Base.Test
using BenchmarkTools
using Grid

tests = ["Fdtd_accuracy", "Fdtd_backprop", "Poisson_testscript_RandomEigenfns", 
	"Models", "Data", ]

println("Running tests:")

for t in tests
	fp = string(t, ".jl")
	println(" * $(fp)")
	include(fp)
end