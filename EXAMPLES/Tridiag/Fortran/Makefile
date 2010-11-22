FC  = gfortran
LD  = gfortran

FORTRAN_FLAGS   = -Wall -g -funderscoring
LD_FLAGS  = -D_THREAD_SAFE=1 -pthread

INCPATH = ../../../INCLUDE
LIBPATH = ../../../LIB

# Libraries to link
LIBS = mrrr gfortran pthread rt m


######################## do not edit below ###########################


FFLAGS   = $(FORTRAN_FLAGS) $(MPIFLAGS) -I$(INCPATH)
LDFLAGS  = $(LD_FLAGS) $(MPIFLAGS) -I$(INCPATH)

.PHONY: all

all: main_all.x main_ind.x main_val.x

# All eigenpairs
main_all.x: maina.o
	$(LD) $(LDFLAGS) $< -L$(LIBPATH) \
        $(foreach LIBRARY,$(LIBS),-l$(LIBRARY)) -o $@

maina.o: maina.f

# Subset of eigenpairs by index
main_ind.x: maini.o
	$(LD) $(LDFLAGS) $< -L$(LIBPATH) \
        $(foreach LIBRARY,$(LIBS),-l$(LIBRARY)) -o $@

maini.o: maini.f

# Subset of eigenpairs by value
main_val.x: mainv.o
	$(LD) $(LDFLAGS) $< -L$(LIBPATH) \
        $(foreach LIBRARY,$(LIBS),-l$(LIBRARY)) -o $@

mainv.o: mainv.f

.PHONY: clean
clean:
	rm -f main_*.x result_*.m core.* *__genmod.* *.o *~
