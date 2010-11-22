FC  = gfortran
LD  = gfortran

FFLAGS   = -Wall -g -funderscoring
LDFLAGS  = -D_THREAD_SAFE=1 -pthread

INCPATH = ../../../INCLUDE
LIBPATH = ../../../LIB

MRRRPATH = ../../../LIB
LAPACK = $(HOME)/libs/lapack-3.2.2/lapack_gnu_LINUX.a
BLAS = $(HOME)/libs/lapack-3.2.2/blas_gnu_LINUX.a

FFLAGS  += -I$(INCPATH)
LDFLAGS += -I$(INCPATH)

main_hermitian.x: main_hermitian.o
	$(LD) $(LDFLAGS) main_hermitian.o -L$(MRRRPATH) -lmrrr \
        $(LAPACK) $(BLAS) -lgfortran -lm -o $@

main_hermitian.o:


.PHONY: clean
clean:
	rm -f main_hermitian.o main_hermitian.x *~ core.*