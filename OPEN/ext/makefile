# Makefile to build extensions

# --- macros
CC = gcc
CFLAGS= -O3

FCC = gfortran
FFLAGS= -ffree-form -O3

FPY = f2py
FPYFLAGS = --f77flags="-O3 -Ofast" --quiet

OBJECTS = hello.o complex-math.o
LIBS = -L/usr/local/lib/sprng/lib -llcg -L/usr/local/lib/pgplot -lcpgplot -lpgplot -lX11 -lftn -lm


# --- targets
all: blombscargle.so 
	@echo "\n Done!"


## compile bayesian lomb scargle periodogram
blombscargle.so: lombscargle_bayes.f90
	@$(FPY) -m blombscargle -c $< $(FPYFLAGS)

# --- remove binary and executable files
clean:
	rm -f blombscargle.so $(OBJECTS)