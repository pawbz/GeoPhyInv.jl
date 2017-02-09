using SeismicInversion

println("testing Julia wrappers for f90 libraries...")

# simple print
ccall( (:sprint, SeismicInversion.f90libs.test), Void, () )

# simple print using OMP
ccall( (:sprint_core, SeismicInversion.f90libs.test), Void, () )


# simple matrix in out
A = rand(4);
B = rand(4);
println(typeof(A))

println("in julia A is ", A)
ccall( (:matrix_in_out, SeismicInversion.f90libs.test), Void, (Ptr{Array{Float64,1}}, Ptr{Float64}), A, B )
println("in julia B is ", B)

# derived data types
#type model


#end #type
