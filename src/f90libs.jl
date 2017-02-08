
module f90libs

# maro to return the absolute path of the library
macro libname(lib)
	return :( joinpath(Pkg.dir("SeismicInversion"), $lib, ".so" ))
end # macro


export libname

end # module 

