# compiler
FC      := gfortran

# compiler flags
FCFLAGS = -ffree-form
# FCFLAGS += -pg -fbacktrace -fno-align-commons 
FCFLAGS += -fopenmp
FCFLAGS += -O3 -Ofast

FPYFLAGS += --opt="-O3" --fcompiler=gnu95 -lgomp
FPYFLAGS += --f77flags="-fopenmp"
FPYFLAGS += --f90flags="-fopenmp"


LIBFLAGS = -L. -lpdflib -lrnglib
# pgplot
PGPLOTFLAGS = -L/opt/mesasdk/pgplot -lpgplot -lX11 -lz
# lapack / blas
LAPACK = -L/usr/lib -llapack -lblas

all: problem1.f90 dream.f90 pgmcmc.f90 get_rvN.f90 pdflib.f90 rnglib.f90 dream_wrapper.f90
	$(FC) -c $(FCFLAGS) dream_wrapper.f90 $(PGPLOTFLAGS) $(LAPACK)
	$(FC) -c $(FCFLAGS) get_rvN.f90 $(PGPLOTFLAGS) $(LAPACK)
	$(FC) -c $(FCFLAGS) problem1.f90 $(PGPLOTFLAGS) $(LAPACK)
	$(FC) -c $(FCFLAGS) dream.f90 $(PGPLOTFLAGS) $(LAPACK)
	$(FC) -c $(FCFLAGS) pgmcmc.f90 $(PGPLOTFLAGS) $(LAPACK)
	$(FC) -c $(FCFLAGS) pdflib.f90 $(PGPLOTFLAGS) $(LAPACK)
	$(FC) -c $(FCFLAGS) rnglib.f90 $(PGPLOTFLAGS) $(LAPACK)
	$(FC) $(FCFLAGS) dream_wrapper.o dream.o problem1.o get_rvN.o pgmcmc.o pdflib.o rnglib.o $(PGPLOTFLAGS) $(LAPACK)


dream.so: problem1.f90 dream.f90 get_rvN.f90 pdflib.f90 rnglib.f90
	f2py -m dream -c $^ $(PGPLOTFLAGS) $(LAPACK) $(FPYFLAGS)

pdflib: pdflib.f90
	$(FC) -c $(FCFLAGS) $^
	ar qc libpdflib.a pdflib.o
	@rm pdflib.o

rnglib: rnglib.f90
	$(FC) -c $(FCFLAGS) $^
	ar qc librnglib.a rnglib.o
	@rm rnglib.o

libs: pdflib rnglib

mcmc: pgmcmc.f90
	$(FC) -o visu_plot $(FCFLAGS) pgmcmc.f90 $(PGPLOTFLAGS)


