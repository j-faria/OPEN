# for compiling python extensions
f2py_version := $(shell f2py -v 2>/dev/null)
ifeq ($(f2py_version), 2)
	FPY = f2py
else
	FPY = f2py2.7
endif

F2PYFLAGS = --quiet
F2PYFLAGS += --opt="-O3"
F2PYFLAGS += --f90flags="-w -Wno-unused-parameter -fPIC -ffree-line-length-none"

# for compiling Fortran
FC := gfortran

FFLAGS += -O3 -ffast-math -fPIC -ffree-line-length-none #-fcheck=all
ifeq ($(PROFFILING),"true")
	FFLAGS += -pg
endif

LIBS = -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB) # $(FLIB)
# LIBS += -L/home/jfaria/data/ranlib.f -lranlib -lblas # this has to go!
LIBS += -L/home/joao/Work/Python/f_ext/ranlib.f -lranlib -lblas # this has to go!

ifndef WITHOUT_MPI
# compiling with MPI support by default. You can disable this with 
# $ make WITHOUT_MPI=1
# but MultiNest must have been compiled without MPI also
FC = mpif90 -f90=${FCPATH}
FFLAGS += -DMPI
endif

MN_OBJ = stringmod.o array_utils.o gp_utils.o
MN_OBJ += params.o like.o nestwrap.o main.o
MN_OBJ += likelihood.o get_rvN.o
# MN_OBJ += 
# MN_OBJ += plotting.o


all:
	echo $(FPY)

# all: nest get_rv1.so get_rvN.so # covmat
fortran: nest

get_rvN.f90: ../ext/get_rvN.f90
	ln -s -f $< $@

test: get_rv1.so get_rvN.so
	python test.py
	ndiff test/out_test.txt test/out_normal.txt -relerr 1.0e-5

get_rv1.so: ../ext/get_rv.f90
	$(FPY) -m get_rv1 -c $^ $(F2PYFLAGS)

get_rvN.so: ../ext/get_rvN.f90 likelihood.f90
	$(FPY) -m get_rvN -c $^ $(F2PYFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) $(LIBS) -I$(NESTLIBDIR) -c -o $@ $^ 
 
nest: $(MN_OBJ)
	$(FC) -o $@ $(MN_OBJ) $(FFLAGS) $(LIBS)


main: main.f90 like.f90 nestwrap.f90 params.f90
	$(FPY) -m main -c $^ $(F2PYFLAGS) $(LIBS)

#covmat: covmat_dump.f90
#	$(FC) -o covmat $(FFLAGS) $^ -L$(FLIBDIR) -I$(FLIBDIR) $(LAPACKLIB) -lmodules

gp.so: gp_python_wrapper.f90 gp_utils.o array_utils.o get_rvN.o
	f2py -c -m gp $^ -llapack $(F2PYFLAGS) -L/home/joao/Work/Python/f_ext/ranlib.f -lranlib -lblas --f90exec="${FCPATH}"
	cp gp.so ../ext/



clean: 
	rm -f nest *.so *.o $(MN_OBJ) *.mod

cleanall:
	rm -f nest *.so *.o *.mod $(MN_OBJ) chains/* covmat

cleanchains:
	rm -f chains/*
