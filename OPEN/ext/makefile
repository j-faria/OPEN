# Makefile to build extensions

# --- macros
CC = gcc
CFLAGS= -O3

FCC = gfortran
FFLAGS= -ffree-form -O3

FPY = f2py
FPYFLAGS = --f77flags="-O3 -Ofast" --quiet
FPYFLAGS += --f90exec="/opt/mesasdk/bin/gfortran" --opt="-O3 -Ofast"
FPYFLAGS += --f90flags="-w -Wno-unused-parameter -fPIC -ffree-line-length-none"

# --- targets
all: glombscargle.so blombscargle.so get_rv.so
	@echo "\n Done!"


## compile bayesian lomb scargle periodogram
blombscargle.so: lombscargle_bayes.f90
	@$(FPY) -m blombscargle -c $< $(FPYFLAGS)

## compile generalized lomb scargle periodogram
glombscargle.so: lombscargle_generalized.f90
	@$(FPY) -m glombscargle -c $< $(FPYFLAGS)

## compile radial velocity curve builder
get_rv.so: get_rv.f90
	@$(FPY) -m get_rv -c $< $(FPYFLAGS)


## compile Dawson-Fabrycky special periodogram
periodogram_DF.so: periodogram_DF.c
	@python setup_periodogram_DF.py build_ext --inplace

# --- remove binary and executable files
clean:
	rm -f blombscargle.so get_rv.so