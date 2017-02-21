
###------------------ some folders -------------------------------
folder=$(shell pwd)/
doc_folder=$(shell pwd)/docs/
bin_folder=$(shell pwd)/bin/
lib_folder=$(shell pwd)/lib/
f90src_folders=$(shell find $(shell pwd) -type f -name "*.f90" | sed -r 's|/[^/]+$$||' |sort |uniq)

include makefile.in
all:
	$(foreach var,$(f90src_folders), make lib_folder=$(lib_folder) -I $(folder) -C $(var) -f $(folder)makef90.mk)
	

.PHONY: docs
docs: 
	julia --color=yes $(doc_folder)/make.jl



# clean *.o, *.d, *.mod  and executable files
.PHONY: clean
clean:
	$(foreach var,$(f90src_folders), rm -f $(var)/*.mod)
	$(foreach var,$(f90src_folders), rm -f $(var)/*.d)
	$(foreach var,$(f90src_folders), rm -f $(var)/*.so)
	$(foreach var,$(f90src_folders), rm -f $(var)/*.o)
	rm -f $(bin_folder)*
	rm -f $(lib_folder)*.so
	rm -f $(lib_folder)*.o

# Show help.
help:
	@echo 'Generic Makefile version 0.0'
	@echo 'Copyright (C) 2016 pawbz@mit.edu'
	@echo
	@echo 'Usage: make [TARGET]'
	@echo 'TARGETS:'
	@echo '  all       compile and link.'
	@echo '  clean     clean.'
	@echo '  help      print help.'
	@echo
	@echo 'Report bugs to pawbz@mit.edu'

# Show variables (for debug use only.)
show:
	@echo 'program                    :' $(program)
	@echo 'f90 source folder(s)       :' $(f90src_folders)
	@echo 'shared lib folder          :' $(lib_folder)
	@echo 'FORTRAN90 complier         :' $(f90comp)
	@echo 'FORTRAN90 compiler flags   :' $(HEADERS)
	@echo 'SOURCES     :' $(SOURCES)
	@echo 'SRC_CXX     :' $(SRC_CXX)
	@echo 'OBJS        :' $(OBJS)
	@echo 'DEPS        :' $(DEPS)
	@echo 'DEPEND      :' $(DEPEND)
	@echo 'COMPILE.c   :' $(COMPILE.c)
	@echo 'COMPILE.cxx :' $(COMPILE.cxx)
	@echo 'LINK.c      :' $(LINK.c)
	@echo 'LINK.cxx    :' $(LINK.cxx)

## End of the makefile
