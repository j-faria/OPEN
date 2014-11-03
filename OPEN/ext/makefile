# Makefile to build extensions

# --- macros
CC = gcc
CFLAGS= -O3

FCC = $(FCPATH)
FFLAGS= -ffree-form -O3

f2py_version := $(shell f2py -v 2>/dev/null)
ifeq ($(f2py_version), 2)
	FPY = f2py
else
	FPY = f2py2.7
endif
FPYFLAGS = --quiet
# FPYFLAGS = --verbose
FPYFLAGS += --opt="-O3" --fcompiler=gnu95 -lgomp
# FPYFLAGS += --f77flags="-fPIC -ffree-line-length-none -fopenmp"
FPYFLAGS += --f90flags="-fPIC -ffree-line-length-none -fopenmp"


SOS = glombscargle.so \
      blombscargle.so \
      get_rv.so \
      get_rvN.so \
      periodogram_DF.so \
      periodogram_CLEAN.so \
      lopast.so \
      julian.so

# --- targets
all: $(SOS)
	@echo "Done!"


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

get_rvN_MultiSite.so: get_rvN_MultiSite.f90
	@$(FPY) -m get_rvN_MultiSite -c $< $(FPYFLAGS)


## compile gp helper functions
GPPATH = ../multinest
gp.so: $(GPPATH)/gp.so
	cp $(GPPATH)/gp.so .

## compile low-pass filter function
lopast.so: lopast.f90
	@$(FPY) -m lopast -c $< $(FPYFLAGS)

## compile Julian date conversion functions
julian.so: julday.f caldat.f
	@$(FPY) -m julian -c $^ $(FPYFLAGS)

## compile Dawson-Fabrycky special periodogram
periodogram_DF.so: periodogram_DF.c
	# echo $(PYTHON)
	@$(PYTHON) setup_periodogram_DF.py build_ext --inplace

## compile CLEAN algorithm periodogram
periodogram_CLEAN.so: periodogram_CLEAN.c
	@$(PYTHON) setup_periodogram_CLEAN.py build_ext --inplace	

# --- remove binary and executable files
clean:
	rm -f $(SOS)


## compile bayesian lomb scargle periodogram with additional outputs
blombscargle_stub.so: lombscargle_bayes_stub.f90
	@$(FPY) -m blombscargle_stub -c $< $(FPYFLAGS)