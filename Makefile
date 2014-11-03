#!/bin/sh

### choose your python interpreter here
PYTHON=/home/joao/anaconda/bin/python
# PYTHON=python2.7

### MultiNest directory (no need to change this)
### relative to OPEN/OPEN/multinest
NESTLIBDIR=../../MultiNest/MultiNest_v3.7

### choose the path to gfortran here
# FCPATH = gfortran
FCPATH = /opt/mesasdk/bin/gfortran

### LAPACK library
LAPACKLIB = -llapack 
# LAPACKLIB = -L/data/jfaria/mesasdk/mesasdk -llapack 

PROFFILING = "false"


NO_COLOR=\033[0m
OK_COLOR=\033[0;32m
ERROR_COLOR=\033[0;31m
WARN_COLOR=\033[0;33m

OK_STRING=$(OK_COLOR)[OK]$(NO_COLOR)
ERROR_STRING=$(ERROR_COLOR)[ERRORS]$(NO_COLOR)
WARN_STRING=$(WARN_COLOR)[WARNING]$(NO_COLOR)


mpif90_version := $(shell mpif90 -v 2>/dev/null)
no_mpif90 = "mpif90 is not installed. MultiNest will be compiled without MPI support"

### export all variables to nested makes
export

all: ext nest

ext: 
	@echo "Compiling OPEN - Fortran extensions..."
	@make -C ./OPEN/ext >> compile.out
	@echo "OPEN extensions  --  $(OK_STRING)"
	@echo

.PHONY: MultiNest
MultiNest: 
	@echo "Compiling MultiNest..."
ifdef mpif90_version
	@echo "$(OK_COLOR) Found $(mpif90_version)$(NO_COLOR)"
	@make -C ./MultiNest/MultiNest_v3.7
else
	@echo "$(WARN_COLOR) $(no_mpif90)$(NO_COLOR)"
	@make -C ./MultiNest/MultiNest_v3.7 WITHOUT_MPI=1
endif
	@make libnest3.so -C ./MultiNest/MultiNest_v3.7 --quiet 
	@echo "MultiNest  --  $(OK_STRING)"
	@echo 
#	might need this
#	sudo cp ./MultiNest/MultiNest_v3.7/libnest3.so /usr/local/lib/

nest: MultiNest
	@echo "Compiling OPEN-MultiNest interface..."
	@make clean -C ./OPEN/multinest --quiet 
ifdef mpif90_version
	make nest -C ./OPEN/multinest
	@make gp.so -C ./OPEN/multinest --quiet
else
	@make nest -C ./OPEN/multinest WITHOUT_MPI=1 --quiet
endif
	@echo "OPEN <-> MultiNest  --  $(OK_STRING)"


gp: 
	@make gp.so -C ./OPEN/multinest

clean: clean-ext clean-multinest clean-multinest-open

clean-ext:
	make clean -C ./OPEN/ext
clean-multinest:
	make clean -C ./OPEN/multinest
clean-multinest-open:
	make clean -C ./MultiNest/MultiNest_v3.7

