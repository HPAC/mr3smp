CC = gcc
LD = gcc

INCPATH = ../../../INCLUDE

CFLAGS  = -Wall -g -I$(INCPATH) -DCOMPLEX_SUPPORTED
LDFLAGS = -Wall -g -pthread -I$(INCPATH)

MRRRPATH = ../../../LIB
LAPACKPATH = $(HOME)/libs/lapack-3.2.2
BLASPATH = $(HOME)/libs/lapack-3.2.2

main_hermitian.x: main_hermitian.o
	$(LD) $(LDFLAGS) main_hermitian.o -L$(MRRRPATH) -L$(LAPACKPATH) -L$(BLASPATH) -lmrrr -llapack_gnu_LINUX -lblas_gnu_LINUX -lgfortran -lm -o $@

main_hermitian.o:


.PHONY: clean
clean:
	rm -f main_hermitian.o main_hermitian.x *~ core.*