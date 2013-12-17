# for compiling python extensions
FPY = f2py
F2PYFLAGS = --quiet
F2PYFLAGS += --opt="-O3 -Ofast"
F2PYFLAGS += --f90flags="-w -Wno-unused-parameter -fPIC -ffree-line-length-none"

# for compiling Fortran
FC := gfortran
FCPATH = /opt/mesasdk/bin/gfortran

FFLAGS +=  -I. -O3 -w -Wno-unused-parameter -fPIC -ffree-line-length-none #-fcheck=all

# utils library directory (CHANGE THIS!!!)
FLIBDIR = /home/joao/Software/fortranlib

# for MultiNest (CHANGE THIS!!!)
NESTLIBDIR=/home/joao/Software/MultiNest_v3.4


FLIB = -L$(FLIBDIR) -I$(FLIBDIR) -lmodules
LAPACKLIB = -llapack -L/usr/lib
LIBS = -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB) $(FLIB)

ifndef WITHOUT_MPI
# compiling with MPI support by default. You can disable this with 
# $ make WITHOUT_MPI=1
# but MultiNest must have been compiled without MPI also
FC = mpif90 -f90=${FCPATH}
FFLAGS += -DMPI
endif

MN_OBJ = params.o like.o nestwrap.o main.o
MN_OBJ += likelihood.o get_rvN.o
# MN_OBJ += plotting.o


all: nest get_rv1.so get_rvN.so # covmat
fortran: nest # covmat

test: get_rv1.so get_rvN.so
	python test.py
	ndiff test/out_test.txt test/out_normal.txt -relerr 1.0e-5

get_rv1.so: get_rv1.f90
	$(FPY) -m get_rv1 -c $^ $(F2PYFLAGS)

get_rvN.so: get_rvN.f90 likelihood.f90
	$(FPY) -m get_rvN -c $^ $(F2PYFLAGS)

%.o %.mod: %.f90
	$(FC) $(FFLAGS) -I$(NESTLIBDIR)  $(LIBS) -c $*.f90
 
nest: $(MN_OBJ)
	echo $(FLIB)
	$(FC) -o nest $(MN_OBJ) $(FFLAGS) $(LIBS)


main: main.f90 like.f90 nestwrap.f90 params.f90
	$(FPY) -m main -c $^ $(F2PYFLAGS) $(LIBS)

covmat: covmat_dump.f90
	$(FC) -o covmat $(FFLAGS) $^ -L$(FLIBDIR) -I$(FLIBDIR) $(LAPACKLIB) -lmodules





clean: 
	rm -f nest *.so *.o *.mod $(MN_OBJ) 

cleanall:
	rm -f nest *.so *.o *.mod $(MN_OBJ) chains/* covmat

cleanchains:
	rm -f chains/*
