#!/bin/sh

### choose your python interpreter here
# PYTHON=/home/joao/anaconda/bin/python
PYTHON=python

### utils library directory (CHANGE THIS!!!)
FLIBDIR = /home/joao/Software/fortranlib

### MultiNest directory 
### relative to OPEN/OPEN/multinest
NESTLIBDIR=../../MultiNest/MultiNest_v3.7

### gfortran path
# FCPATH = /opt/mesasdk/bin/gfortran
FCPATH = /data/jfaria/mesasdk/mesasdk/bin/gfortran

### LAPACK library
LAPACKLIB = -L/data/jfaria/mesasdk/mesasdk/lib64 -llapack 
# LAPACKLIB = -L/data/jfaria/mesasdk/mesasdk -llapack 


NO_COLOR=\033[0m
OK_COLOR=\033[92m
ERROR_COLOR=\033[91m
WARN_COLOR=\033[93m

OK_STRING=$(OK_COLOR)[OK]$(NO_COLOR)
ERROR_STRING=$(ERROR_COLOR)[ERRORS]$(NO_COLOR)
WARN_STRING=$(WARN_COLOR)[WARNING]$(NO_COLOR)


mpif90_version := $(shell mpif90 -v 2>/dev/null)
no_mpif90 = "mpif90 is not installed. MultiNest will be compiled without MPI support"

### export all variables to nested makes
export

all: ext nest

ext: 
	@echo "Compiling Fortran extensions..."
	make -C ./OPEN/ext
	@echo "OPEN extensions  --  $(OK_STRING)"


.PHONY: MultiNest
MultiNest: 
	@echo "Compiling MultiNest..."
#	# @command -v mpif90 >/dev/null 2>&1 || { echo >&2 "$(WARN_STRING) $(no_mpif90)"; }
ifdef mpif90_version
	@echo "$(OK_COLOR) Found $(mpif90_version)$(NO_COLOR)"
	@make -C ./MultiNest/MultiNest_v3.7
else
	@echo "$(WARN_COLOR) $(no_mpif90)$(NO_COLOR)"
	@make -C ./MultiNest/MultiNest_v3.7 WITHOUT_MPI=1
endif
	@make libnest3.so -C ./MultiNest/MultiNest_v3.7 --quiet 
	@echo "MultiNest  --  $(OK_STRING)"


nest: MultiNest
	@echo "Compiling OPEN-MultiNest interface..."
	@make clean -C ./OPEN/multinest --quiet 
ifdef mpif90_version
	@make nest -C ./OPEN/multinest
else
	@make nest -C ./OPEN/multinest WITHOUT_MPI=1
endif
	@echo "OPEN <-> MultiNest  --  $(OK_STRING)"


clean: clean-ext clean-multinest clean-multinest-open

clean-ext:
	make clean -C ./OPEN/ext
clean-multinest:
	make clean -C ./OPEN/multinest
clean-multinest-open:
	make clean -C ./MultiNest/MultiNest_v3.7
