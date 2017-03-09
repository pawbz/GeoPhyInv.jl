
module F90libs

# macro to return the absolute path of the library
function libname(lib)
    libname = string(lib , ".so");
    return ( joinpath(Pkg.dir("SeismicInversion"), "lib", libname))
end # macro

macro make_global_const(lib, libname)
	return :(global const $lib = string($libname))
end

# function to make global constants of library names
function __init__()
#	for lib in readdir(joinpath(Pkg.dir("SeismicInversion"), "lib"))
#		println("importing f90 library", lib)
#		@make_global_const(lib, string(lib))
#	end

global const test = libname("test") 
global const fdtd = libname("fdtd") 
global const interpolation = libname("interpolation") 
global const io = libname("io") 
global const str = libname("string_routines") 
global const fft = libname("fft") 


end


end # module 

