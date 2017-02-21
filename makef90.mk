# *.mod are module files
# *.d are the files containing depedencies for each module file
# the *.d files in C can normally be generated using the compiler option
# '-M'. But I cannot find such an option for gfrotran .. maybe it exists 
# for latest ifort compiler.
# Here I try to prepare *.d files, by trying to grep all the lines containing
# use statement. 
# THe modules should contain one of the following two statements:
#
#  1) use module_filename, only : foo! USEINDEPFILE
#  2) use module_filename! USEINDEPFILE
#
#  3) use module_name
#     !INCLUDEINDEPFILE module_filename,
# 
# If the module name is same as the filename containing the module, u can 
# use the first two statements. Otherwise the third method is a general case
# where a dependency on 'module_filename' is added in the .d file.
#
#
#
# F90 makefile created by Pawan Bharadwaj, 13 July'14
# 10 Jan'15 included rtabs to remove tab spaces from f90 files

# All the files in $PWD with .f90 extension are compiled. 
# If you want any file not to be compiled, please use a different extension,
# for example .F90.
#
#
include makefile.in
# Names of object, dependency files:
f90files = $(shell find . -name "*.f90")
all_objs = $(f90files:./%.f90=%.o)
all_shared_libs = $(f90files:./%.f90=%.so)
all_deps = $(f90files:./%.f90=%.d)
remove_tabs = $(f90files:./%.f90=%.rtabs)


all: $(all_shared_libs)
	@echo $(program)
	@echo 'shared libraries of f90 files generated successfully'

show: 
	@echo 'program                    :' $(f90files)
	@echo 'program                    :' $(all_objs)
	@echo 'program                    :' $(all_shared_libs)
	@echo 'program                    :' $(all_deps)
	@echo 'program                    :' $(remove_tabs)
	

# include *.d files that are generated
include $(all_deps)

# generate *.so, shared libraries for all the modules
%.so: %.f90 #%.rtabs
	$(f90comp) $(openmp) $(debug) $(gprof) $(speedup) $(static) $(double) $(fftlib) $(blas) $(saclib) \
		-fPIC -shared $(shell echo $^ | sed 's/.so/.f90/g') -o $(lib_folder)$@

# remove tab spaces in fortran source files
rtabs: $(remove_tabs)
%.rtabs: %.f90
	cp $< $*.d
	expand $*.d >$< 
	rm *.d

# generate only *.d files
deps: $(all_deps)
%.d: %.f90
	grep -h 'USEINDEPFILE\|INCLUDEINDEPFILE' $< | \
	sed 's/!/,/' | \
	sed 's/,INCLUDEINDEPFILE/use/' |\
	sed 's/,.*/.so/' | \
	sed 's/use /$*.so: /' \
	>$*.d

## End of the f90 makefile
