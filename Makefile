#!/bin/sh

NO_COLOR=\033[0m
OK_COLOR=\033[92m
ERROR_COLOR=\033[91m
WARN_COLOR=\033[93m

OK_STRING=$(OK_COLOR)[OK]$(NO_COLOR)
ERROR_STRING=$(ERROR_COLOR)[ERRORS]$(NO_COLOR)
WARN_STRING=$(WARN_COLOR)[WARNING]$(NO_COLOR)

mpif90_version := $(shell mpif90 -v 2>/dev/null)
no_mpif90 = "mpif90 is not installed. MultiNest will be compiled without MPI support"

all: ext nest


ext: 
	@echo "Compiling Fortran extensions..."
	@make -C ./OPEN/ext >> compile.out
	@echo "OPEN extensions  --  $(OK_STRING)"




.PHONY: MultiNest
MultiNest: 
	@echo "Compiling MultiNest..."
#	# @command -v mpif90 >/dev/null 2>&1 || { echo >&2 "$(WARN_STRING) $(no_mpif90)"; }
ifdef mpif90_version
	@echo "$(OK_COLOR) Found $(mpif90_version)$(NO_COLOR)"
	@make -C ./MultiNest/MultiNest_v3.7 --quiet 
else
	@echo "$(WARN_COLOR) $(no_mpif90)$(NO_COLOR)"
	@make -C ./MultiNest/MultiNest_v3.7 WITHOUT_MPI=1 --quiet
endif
	@make libnest3.so -C ./MultiNest/MultiNest_v3.7 --quiet 
	@echo "MultiNest  --  $(OK_STRING)"


nest: MultiNest
	@echo "Compiling OPEN-MultiNest interface..."
	@make clean -C ./OPEN/multinest --quiet 
ifdef mpif90_version
	@make fortran -C ./OPEN/multinest --quiet 
else
	@make fortran -C ./OPEN/multinest WITHOUT_MPI=1 --quiet 
	@echo "OPEN <-> MultiNest  --  $(OK_STRING)"


clean:
	make clean -C ./OPEN/ext
	make clean -C ./OPEN/multinest
	make clean -C ./MultiNest/MultiNest_v3.7
