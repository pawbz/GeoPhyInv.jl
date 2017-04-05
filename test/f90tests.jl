using SIT
using Base.Test

println("testing Julia wrappers for f90 libraries...")

# simple print
println(repeat(" ",5), "simple print...")
sayhelloin = "Julia"
sayhelloout = "Jua"
ccall( (:sprint, SIT.F90libs.test), Void, (Ptr{UInt8},Ptr{UInt8}), sayhelloin, sayhelloout )
println(sayhelloout)
#@test isequal(sayhelloin, sayhelloout)

# simple print using OMP
println(repeat(" ",5), "simple print using OMP...")
ccall( (:sprint_core, SIT.F90libs.test), Void, () )


# simple matrix in out
println(repeat(" ",5), "simple vectors in and out..")
n1in = 4; n2in = 2;
Ain = rand(n1in, n2in);
A = rand(n1in, n2in);
kin = [0.001]; 
k = [0.02]; n1 = [0]; 

ccall( (:matrix_in_out, SIT.F90libs.test), Void, (Ptr{Float64}, Ptr{Float64}, Ref{Int64}, 
							       Ref{Int64}), Ain, A, n1in, n2in)
@test_approx_eq A Ain 

ccall( (:real64_in_out, SIT.F90libs.test), Void, (Ptr{Float64},Ptr{Float64}), kin, k)
@test_approx_eq k kin 

ccall( (:int64_in_out, SIT.F90libs.test), Void, (Ref{Int64},Ref{Int64}), n1in, n1)
n1 = n1[1]
@test_approx_eq n1 n1in 

