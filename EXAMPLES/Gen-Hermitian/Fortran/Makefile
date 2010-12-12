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

all: main_gherm_all.x main_gherm_ind.x

main_gherm_all.x: main_gherm_all.o
	$(LD) $(LDFLAGS) main_gherm_all.o -L$(MRRRPATH) -lmrrr \
        $(LAPACK) $(BLAS) -lgfortran -lm -o $@

main_gherm_all.o:

main_gherm_ind.x: main_gherm_ind.o
	$(LD) $(LDFLAGS) main_gherm_ind.o -L$(MRRRPATH) -lmrrr \
        $(LAPACK) $(BLAS) -lgfortran -lm -o $@

main_gherm_ind.o:

.PHONY: clean
clean:
	rm -f main_*.o main_*.x *~ core.*