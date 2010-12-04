FC  = gfortran
LD  = gfortran

FFLAGS   = -Wall -g -funderscoring
LDFLAGS  = -D_THREAD_SAFE=1 -pthread

INCPATH = ../../../INCLUDE
LIBPATH = ../../../LIB 

# # Build without including LAPACK and BLAS routines  
# # by setting INCLAPACK to zero
# LAPACKPATH = $(HOME)/libs/lapack-3.3.0
# LIBS = mrrr lapack_gnu_LINUX blas_gnu_LINUX gfortran m pthread rt

# Otherwise
LAPACKPATH = 
LIBS = mrrr gfortran m pthread rt


######################## do not edit below ###########################


FFLAGS  += -I$(INCPATH)
LDFLAGS += -I$(INCPATH)
LIBPATH += $(LAPACKPATH)


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
