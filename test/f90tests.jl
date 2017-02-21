using SeismicInversion

println("testing Julia wrappers for f90 libraries...")

# simple print
println(repeat(" ",5), "simple print...")
sayhello = "Julia"
ccall( (:sprint, SeismicInversion.f90libs.test), Void, (Ptr{UInt8},), sayhello )

# simple print using OMP
println(repeat(" ",5), "simple print using OMP...")
ccall( (:sprint_core, SeismicInversion.f90libs.test), Void, () )


# simple matrix in out
println(repeat(" ",5), "simple vectors in and out..")
n1 = 4; n2 = 2;
A = rand(n1, n2);
B = rand(n1);

println(repeat(" ",10), "in julia A is ", A)
println(repeat(" ",10), "in julia input B is ", B)
ccall( (:matrix_in_out, SeismicInversion.f90libs.test), Void, (Ptr{Float64}, Ptr{Float64}, Ref{Int64}, 
							       Ref{Int64}), A, B, n1, n2)
println(repeat(" ",10), "in julia output B is ", B)



# derived data types
#=
n1=4; n2=2;
type Model
	A::Vector{Float64}
end #type
model_in = Model([3,4,5])
ccall( (:matrix_in_out, SeismicInversion.f90libs.test), Void, (Ptr{Float64}, Ptr{Float64}, Ref{Int64}, 
							       Ref{Int64}), A, B, n1, n2)
							       =#
							       
