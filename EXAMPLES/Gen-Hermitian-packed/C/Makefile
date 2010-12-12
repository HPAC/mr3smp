CC = gcc
LD = gcc

INCPATH = ../../../INCLUDE

CFLAGS  = -Wall -g -I$(INCPATH) -DCOMPLEX_SUPPORTED
LDFLAGS = -Wall -g -pthread -I$(INCPATH)

MRRRPATH = ../../../LIB
LAPACKPATH = $(HOME)/libs/lapack-3.2.2
BLASPATH = $(HOME)/libs/lapack-3.2.2

all: main_ghp_all.x main_ghp_ind.x

# Computing all eigenpairs
main_ghp_all.x: main_ghp_all.o
	$(LD) $(LDFLAGS) main_ghp_all.o -L$(MRRRPATH) -L$(LAPACKPATH) -L$(BLASPATH) -lmrrr -llapack_gnu_LINUX -lblas_gnu_LINUX -lgfortran -lm -o $@

main_ghp_all.o:

# Computing a subset by index 
main_ghp_ind.x: main_ghp_ind.o
	$(LD) $(LDFLAGS) main_ghp_ind.o -L$(MRRRPATH) -L$(LAPACKPATH) -L$(BLASPATH) -lmrrr -llapack_gnu_LINUX -lblas_gnu_LINUX -lgfortran -lm -o $@

main_ghp_ind.o:


.PHONY: clean
clean:
	rm -f main_ghp_*.o main_ghp_*.x *~ core.*