CC = gcc
LD = gcc

INCPATH = ../../../INCLUDE

CFLAGS  = -Wall -g -I$(INCPATH)
LDFLAGS = -Wall -g -pthread -I$(INCPATH)

MRRRPATH = ../../../LIB
LAPACK = $(HOME)/libs/lapack-3.2.2/lapack_gnu_LINUX.a
BLAS = $(HOME)/libs/lapack-3.2.2/blas_gnu_LINUX.a

main_sym.x: main_sym.o
	$(LD) $(LDFLAGS) main_sym.o -L$(MRRRPATH) -lmrrr $(LAPACK) $(BLAS) \
        -lgfortran -lm -o $@

main_sym.o:


.PHONY: clean
clean:
	rm -f main_sym.o main_sym.x *~ core.*