CC = gcc
LD = gcc

INCPATH = ../../../INCLUDE

CFLAGS  = -Wall -g -I$(INCPATH)
LDFLAGS = -Wall -g -pthread -I$(INCPATH)

MRRRPATH = ../../../LIB
LAPACK = $(HOME)/libs/lapack-3.2.2/lapack_gnu_LINUX.a
BLAS = $(HOME)/libs/lapack-3.2.2/blas_gnu_LINUX.a

main_hermitian.x: main_hermitian.o
	$(LD) $(LDFLAGS) main_hermitian.o -L$(MRRRPATH) -lmrrr \
        $(LAPACK) $(BLAS) -lgfortran -lm -o $@

main_hermitian.o:


.PHONY: clean
clean:
	rm -f main_hermitian.o main_hermitian.x *~ core.*