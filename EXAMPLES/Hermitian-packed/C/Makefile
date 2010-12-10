CC = gcc
LD = gcc

INCPATH = ../../../INCLUDE

CFLAGS  = -Wall -g -I$(INCPATH) -DCOMPLEX_SUPPORTED
LDFLAGS = -Wall -g -pthread -I$(INCPATH)

MRRRPATH = ../../../LIB
LAPACK = $(HOME)/libs/lapack-3.2.2/lapack_gnu_LINUX.a
BLAS = $(HOME)/libs/lapack-3.2.2/blas_gnu_LINUX.a

all: main_all.x main_ind.x

main_all.x: main_hermp_all.o
	$(LD) $(LDFLAGS) main_hermp_all.o -L$(MRRRPATH) -lmrrr $(LAPACK) $(BLAS) -lgfortran -lm -o $@

main_hermp_all.o:

main_ind.x: main_hermp_ind.o
	$(LD) $(LDFLAGS) main_hermp_ind.o -L$(MRRRPATH) -lmrrr $(LAPACK) $(BLAS) -lgfortran -lm -o $@

main_hermp_ind.o:

.PHONY: clean
clean:
	rm -f main_*.o main_*.x *~ core.*