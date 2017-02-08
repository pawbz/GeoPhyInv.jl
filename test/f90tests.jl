using SeismicInversion

println("testing Julia wrappers for f90 libraries...")

f90libname = SeismicInversion.f90libs.@libname("test")

println("simple call...", f90libname)
#ccall( (:F, f90libname, Void, () )
