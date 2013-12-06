# for compiling python extensions
FPY = f2py
F2PYFLAGS = --f90exec="/opt/mesasdk/bin/gfortran" --opt="-O3 -Ofast"
F2PYFLAGS += --f90flags="-w -Wno-unused-parameter -fPIC -ffree-line-length-none"

# for compiling Fortran
FC := gfortran # gfortran by default
FCPATH = /opt/mesasdk/bin/gfortran

FFLAGS +=  -I. -O3 -w -Wno-unused-parameter -fPIC -ffree-line-length-none

# utils library directory
FLIBDIR = /home/joao/Software/fortranlib
FLIB = -L$(FLIBDIR) -I$(FLIBDIR) -lmodules

# for MultiNest
NESTLIBDIR=/home/joao/Software/MultiNest_v3.4
LAPACKLIB = -llapack -L/usr/lib

LIBS = -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB) $(FLIB)

ifndef WITHOUT_MPI
# compiling with MPI support. 
# You can disable this with 
# $ make WITHOUT_MPI=1
FC = mpif90 -f90=${FCPATH}
FFLAGS += -DMPI
endif

MN_OBJ = params.o like.o nestwrap.o main.o
MN_OBJ += likelihood.o get_rvN.o
# MN_OBJ += plotting.o


all: gaussian get_rv1.so get_rvN.so # covmat
fortran: gaussian # covmat

test: get_rv1.so get_rvN.so
	python test.py
	ndiff test/out_test.txt test/out_normal.txt -relerr 1.0e-5

get_rv1.so: get_rv1.f90
	$(FPY) -m get_rv1 -c $^ $(F2PYFLAGS)

get_rvN.so: get_rvN.f90 likelihood.f90
	$(FPY) -m get_rvN -c $^ $(F2PYFLAGS)

# plotting.o plotting.mod: plotting.f90
# 	$(FC) $(FFLAGS) -I$(NESTLIBDIR)  $(LIBS) -c $*.f90

%.o %.mod: %.f90
	$(FC) $(FFLAGS) -I$(NESTLIBDIR)  $(LIBS) -c $*.f90
 
gaussian: $(MN_OBJ)
	echo $(FLIB)
	$(FC) -o gaussian $(MN_OBJ) $(FFLAGS) $(LIBS)


main: main.f90 like.f90 nestwrap.f90 params.f90
	$(FPY) -m main -c $^ $(F2PYFLAGS) $(LIBS)

covmat: covmat_dump.f90
	$(FC) -o covmat $(FFLAGS) $^ -L$(FLIBDIR) -I$(FLIBDIR) $(LAPACKLIB) -lmodules





clean: 
	rm -f gaussian *.so *.o *.mod $(MN_OBJ) 

cleanall:
	rm -f gaussian *.so *.o *.mod $(MN_OBJ) chains/* covmat

cleanchains:
	rm -f chains/*
