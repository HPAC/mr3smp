CC  = gcc
LD  = gcc

CFLAGS   = -Wall -g -pthread
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


CFLAGS  += -I$(INCPATH)
LIBPATH += $(LAPACKPATH)

.PHONY: all

all: main_all.x main_ind.x main_val.x

# All eigenpairs
main_all.x: main_all.o
	$(LD) $(LDFLAGS) $< $(foreach LIBP,$(LIBPATH),-L$(LIBP)) \
        $(foreach LIBRARY,$(LIBS),-l$(LIBRARY)) -o $@

main_all.o: main_all.c

# Subset of eigenpairs by index
main_ind.x: main_ind.o
	$(LD) $(LDFLAGS) $< $(foreach LIBP,$(LIBPATH),-L$(LIBP)) \
        $(foreach LIBRARY,$(LIBS),-l$(LIBRARY)) -o $@

main_ind.o: main_ind.c

# Subset of eigenpairs by value
main_val.x: main_val.o
	$(LD) $(LDFLAGS) $< $(foreach LIBP,$(LIBPATH),-L$(LIBP)) \
        $(foreach LIBRARY,$(LIBS),-l$(LIBRARY)) -o $@

main_val.o: main_val.c

.PHONY: clean
clean:
	rm -f main_*.x core.* *__genmod.* *.o *~
