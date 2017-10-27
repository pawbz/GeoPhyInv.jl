
tests = ["Fdtd_accuracy", "Fdtd_backprop", "Poisson_testscript_RandomEigenfns", "Smooth", "Decon", "Misfits", "DSP",
	"Interpolation", "Models", ]

println("Running tests:")

for t in tests
	fp = string(t, ".jl")
	println(" * $(fp)")
	include(fp)
end
