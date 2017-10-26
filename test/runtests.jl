
tests = ["Fdtd_accuracy", "Poisson_testscript_RandomEigenfns", "Smooth"]

println("Running tests:")

for t in tests
	fp = string(t, ".jl")
	println(" * $(fp)")
	include(fp)
end
