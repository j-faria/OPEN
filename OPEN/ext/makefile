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
all: blombscargle.so get_rv.so
	@echo "\n Done!"


## compile bayesian lomb scargle periodogram
blombscargle.so: lombscargle_bayes.f90
	@$(FPY) -m blombscargle -c $< $(FPYFLAGS)

## compile radial velocity curve builder
get_rv.so: get_rv.f90
	@$(FPY) -m get_rv -c $< $(FPYFLAGS)

# --- remove binary and executable files
clean:
	rm -f blombscargle.so get_rv.so