using GeoPhyInv # precompile it

"""
* `files` are in folder
* `gfiles` are files in the generated folder
"""
function run_test(files, folder, gfiles, nfail)
    exename = joinpath(Sys.BINDIR, Base.julia_exename())
    testdir = @__DIR__
    for file in files
        fp = joinpath(testdir, folder, string(file, ".jl"))
        printstyled("Testing GeoPhyInv file:\t", fp, "\n"; bold = true, color = :white)
        # if f ∈ excludedfiles
        # println("Test Skip:")
        # println("$f")
        # continue
        # end
        try
            run(`$exename -O3 --startup-file=no --check-bounds=no $fp`)
        catch ex
            nfail += 1
        end
    end

    for file in gfiles
        fpg = joinpath(testdir, "generated", folder, string(file, ".jl"))
        printstyled("Testing GeoPhyInv file:\t", fpg, "\n"; bold = true, color = :white)
        try
            run(`$exename -O3 --startup-file=no --check-bounds=no $fpg`)
        catch ex
            nfail += 1
        end
    end
    printstyled("Number of failed tests:\t", nfail, "\n"; bold = true, color = :white)
end


nfail = 0
run_test(["convert", "param_adj"], "media", ["doc"], nfail)
run_test([], "ageom", ["doc"], nfail)
run_test([], "srcwav", ["doc"], nfail)
run_test([], "data", ["doc"], nfail)
# run_test(["accuracy2D", "rho_projection"],"fdtd", ["doc","create_snaps","reuse_expt"], nfail)
#run_test(["accuracy2D","backprop", "rho_projection"],"fdtd", ["doc","create_snaps","reuse_expt"], nfail)
#run_test(["gradient_accuracy", "born_map"],"fwi", ["doc", "pizza", "born_tutorial"], nfail)
# run_test(["interp_tests"], "Interpolation", [], nfail)
# run_test(
# ["testscript_RandomEigenfns", "adj_state_expt", "adj_state", "testdAdx"],
# "Poisson",
# ["doc", "forw", "test_born"],
# nfail,
# )

# exit(nfail)









