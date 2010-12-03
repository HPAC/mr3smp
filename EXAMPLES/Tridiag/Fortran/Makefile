FC  = gfortran
LD  = gfortran

FFLAGS   = -Wall -g -funderscoring
LDFLAGS  = -D_THREAD_SAFE=1 -pthread

# LAPACK and BLAS only needs to be linked if build without including 
# the routines by setting INCLAPACK to zero
INCPATH = ../../../INCLUDE
LIBPATH = ../../../LIB $(HOME)/libs/lapack-3.3
LIBS = mrrr lapack_gnu_LINUX blas_gnu_LINUX gfortran m pthread rt


######################## do not edit below ###########################


FFLAGS  += -I$(INCPATH)
LDFLAGS += -I$(INCPATH)

.PHONY: all

all: main_all.x main_ind.x main_val.x

# All eigenpairs
main_all.x: maina.o
	$(LD) $(LDFLAGS) $< $(foreach LIBP,$(LIBPATH),-L$(LIBP)) \
        $(foreach LIBRARY,$(LIBS),-l$(LIBRARY)) -o $@

maina.o: maina.f

# Subset of eigenpairs by index
main_ind.x: maini.o
	$(LD) $(LDFLAGS) $< $(foreach LIBP,$(LIBPATH),-L$(LIBP)) \
        $(foreach LIBRARY,$(LIBS),-l$(LIBRARY)) -o $@

maini.o: maini.f

# Subset of eigenpairs by value
main_val.x: mainv.o
	$(LD) $(LDFLAGS) $< $(foreach LIBP,$(LIBPATH),-L$(LIBP)) \
        $(foreach LIBRARY,$(LIBS),-l$(LIBRARY)) -o $@

mainv.o: mainv.f

.PHONY: clean
clean:
	rm -f main_*.x result_*.m core.* *__genmod.* *.o *~
