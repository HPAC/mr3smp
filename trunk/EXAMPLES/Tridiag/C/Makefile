CC  = gcc
LD  = gcc

CFLAGS   = -Wall -g -pthread
LDFLAGS  = -D_THREAD_SAFE=1 -pthread

# LAPACK and BLAS only needs to be linked if build without including 
# the routines by setting INCLAPACK to zero
INCPATH = ../../../INCLUDE
LIBPATH = ../../../LIB $(HOME)/libs/lapack-3.3
LIBS = mrrr lapack_gnu_LINUX blas_gnu_LINUX gfortran m pthread rt


######################## do not edit below ###########################


CFLAGS  += -I$(INCPATH)

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
