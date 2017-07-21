using BinDeps

#if !(is_linux())
#	error("currently, SIT only supports Linux")
#end
#
## get the current path
#currentFilePath = @__FILE__()
#currentDirPath = dirname(currentFilePath)
#
#srcf90dir=joinpath(currentDirPath, "src", "f90")
#makefile=joinpath(currentDirPath,"makefile")
#
#
## use mkdir to create usr/lib
#function mklibdir()
#	usrdir = joinpath(currentDirPath, "usr");
#	if isdir(usrdir) == false
#		mkdir(usrdir);
#	end
#
#	libdir = joinpath(usrdir, "lib");
#	if isdir(libdir) == false
#		mkdir(libdir);
#	end
#end
#
#
#f90cmd=string("ifort -assume byterecl  -fp-model fast=2 -Ofast  -autodouble -fPIC -shared")
#f90cmd=string("gfortran -fopenmp   -O3  -fdefault-real-8 -fdefault-double-8 -fPIC -shared")
#
#libsrcvec=[
#	"precision_mod.f90", "string_routines.f90 precision_mod.f90", 
#	"error_messages.f90 precision_mod.f90 string_routines.f90", 
#	"unique.f90 error_messages.f90",
#	"sorting.f90",
#	"io.f90 precision_mod.f90 error_messages.f90 unique.f90 sorting.f90 string_routines.f90",
#	"test.f90 precision_mod.f90 string_routines.f90"
#	]
#
#libsrcvec=[
#	"precision_mod.f90", "string_routines.f90", 
#	"error_messages.f90", 
#	"unique.f90",
#	"sorting.f90",
#	"io.f90",
#	"test.f90"
#	]
#
#libsrcname=broadcast(libsrcvec) do x
#	joinpath(srcf90dir, x)
#end
#
#
#libsovec=[
#	"precision_mod.so",
#	"string_routines.so",
#	"error_messages.so",
#	"unique.so",
#	"sorting.so",
#	"io.so",
#	"test.so"
#	]
#
#libsovec=["lib.so"]
# 
#function compile()
#	for ilib in 1:length(libsovec)
#		libsrcname=joinpath(srcf90dir,libsrcvec[ilib])
#		libsoname = joinpath(currentDirPath, "usr", "lib", libsovec[ilib])
#		run(`gfortran -O3  -fdefault-real-8 -fdefault-double-8 -fPIC -shared $libsrcname -o $libsoname`)
#		println(libsrcname, libsoname)
#
#		outputfile = open(joinpath(currentDirPath, "deps.jl"), "w");
#		write(outputfile, string("@checked_lib ", libsovec[ilib], " ", libsoname))
#		close(outputfile)
#	end
#end
#
#
#function writeDeps()
#	libpath = joinpath(currentDirPath, "usr", "lib", "liblbfgsbf.so")
#	outputfile = open(joinpath(currentDirPath, "deps.jl"), "w");
#	write(outputfile, "macro checked_lib(libname, path)\n")    
#	write(outputfile, "(dlopen_e(path) == C_NULL) && \(error(\"Unable to load \$libname (\$path)")
#	write(outputfile, "Please re-run Pkg.build(package), and restart Julia.\")\n \)")
#        write(outputfile, "quote const \$(esc(libname)) = \$path end\nend\n")
#	close(outputfile)
#end
#
#
##mklibdir()
##compile()
#writeDeps()

