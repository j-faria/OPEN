# Makefile to build extensions

# --- macros
CC = gcc
CFLAGS= -O3

FCC = gfortran
FFLAGS= -ffree-form -O3

FPY = f2py
FPYFLAGS = --quiet
FPYFLAGS += --opt="-O3" --fcompiler=gnu95
# FPYFLAGS += --f77flags="-fPIC -ffree-line-length-none"
FPYFLAGS += --f90flags="-fPIC -ffree-line-length-none"

PYTHON=/home/joao/anaconda/bin/python

# --- targets
all: glombscargle.so blombscargle.so get_rv.so get_rvN.so periodogram_DF.so
	@echo "\n Done!"


## compile bayesian lomb scargle periodogram
blombscargle.so: lombscargle_bayes.f90
	@$(FPY) -m blombscargle -c $< $(FPYFLAGS)

## compile generalized lomb scargle periodogram
glombscargle.so: lombscargle_generalized.f90
	@$(FPY) -m glombscargle -c $< $(FPYFLAGS)

## compile radial velocity curve builders
get_rv.so: get_rv.f90
	@$(FPY) -m get_rv -c $< $(FPYFLAGS)

get_rvN.so: get_rvN.f90
	@$(FPY) -m get_rvN -c $< $(FPYFLAGS)


## compile Dawson-Fabrycky special periodogram
periodogram_DF.so: periodogram_DF.c
	@$(PYTHON) setup_periodogram_DF.py build_ext --inplace

# --- remove binary and executable files
clean:
	rm -f blombscargle.so glombscargle.so get_rv.so get_rvN.so periodogram_DF.so
