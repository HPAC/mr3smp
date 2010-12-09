CC = gcc
LD = gcc

INCPATH = ../../../INCLUDE

CFLAGS  = -Wall -g -I$(INCPATH) -DCOMPLEX_SUPPORTED
LDFLAGS = -Wall -g -pthread -I$(INCPATH)

MRRRPATH = ../../../LIB
LAPACKPATH = $(HOME)/libs/lapack-3.2.2
BLASPATH = $(HOME)/libs/lapack-3.2.2

all: main_genherm_all.x main_genherm_ind.x

# Computing all eigenpairs
main_genherm_all.x: main_genherm_all.o
	$(LD) $(LDFLAGS) main_genherm_all.o -L$(MRRRPATH) -L$(LAPACKPATH) -L$(BLASPATH) -lmrrr -llapack_gnu_LINUX -lblas_gnu_LINUX -lgfortran -lm -o $@

main_genherm_all.o:

# Computing a subset by index 
main_genherm_ind.x: main_genherm_ind.o
	$(LD) $(LDFLAGS) main_genherm_ind.o -L$(MRRRPATH) -L$(LAPACKPATH) -L$(BLASPATH) -lmrrr -llapack_gnu_LINUX -lblas_gnu_LINUX -lgfortran -lm -o $@

main_genherm_ind.o:


.PHONY: clean
clean:
	rm -f main_genherm_*.o main_genherm_*.x *~ core.*