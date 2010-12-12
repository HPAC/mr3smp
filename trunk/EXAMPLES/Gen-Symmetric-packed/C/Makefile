CC = gcc
LD = gcc

INCPATH = ../../../INCLUDE

CFLAGS  = -Wall -g -I$(INCPATH)
LDFLAGS = -Wall -g -pthread -I$(INCPATH)

MRRRPATH = ../../../LIB
LAPACK = $(HOME)/libs/lapack-3.2.2/lapack_gnu_LINUX.a
BLAS = $(HOME)/libs/lapack-3.2.2/blas_gnu_LINUX.a

all: main_gsymp_all.x main_gsymp_ind.x

main_gsymp_all.x: main_gsymp_all.o
	$(LD) $(LDFLAGS) main_gsymp_all.o -L$(MRRRPATH) -lmrrr $(LAPACK) $(BLAS) -lgfortran -lm -o $@

main_gsymp_all.o:

main_gsymp_ind.x: main_gsymp_ind.o
	$(LD) $(LDFLAGS) main_gsymp_ind.o -L$(MRRRPATH) -lmrrr $(LAPACK) $(BLAS) -lgfortran -lm -o $@

main_gsymp_ind.o:

.PHONY: clean
clean:
	rm -f main_gsymp_*.o main_gsymp_*.x *~ core.*