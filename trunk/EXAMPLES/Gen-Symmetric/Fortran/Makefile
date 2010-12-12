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

all: main_gsym_all.x main_gsym_ind.x

main_gsym_all.x: main_gsym_all.o
	$(LD) $(LDFLAGS) main_gsym_all.o -L$(MRRRPATH) -lmrrr $(LAPACK) $(BLAS) -lgfortran -lm -o $@

main_gsym_all.o:

main_gsym_ind.x: main_gsym_ind.o
	$(LD) $(LDFLAGS) main_gsym_ind.o -L$(MRRRPATH) -lmrrr $(LAPACK) $(BLAS) -lgfortran -lm -o $@

main_gsym_ind.o:

.PHONY: clean
clean:
	rm -f main_*.o main_*.x *~ core.*