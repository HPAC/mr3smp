CC  = gcc
LD  = gcc

C_FLAGS   = -Wall -g -pthread
LD_FLAGS  = -D_THREAD_SAFE=1 -pthread

INCPATH = ../../../INCLUDE
LIBPATH = ../../../LIB
LIBS = mrrr gfortran m pthread rt


######################## do not edit below ###########################


CFLAGS   = $(C_FLAGS) -I$(INCPATH)
LDFLAGS  = $(LD_FLAGS)

.PHONY: all

all: main_all.x main_ind.x main_val.x

# All eigenpairs
main_all.x: main_all.o
	$(LD) $(LDFLAGS) $< -L$(LIBPATH) \
        $(foreach LIBRARY,$(LIBS),-l$(LIBRARY)) -o $@

main_all.o: main_all.c

# Subset of eigenpairs by index
main_ind.x: main_ind.o
	$(LD) $(LDFLAGS) $< -L$(LIBPATH) \
        $(foreach LIBRARY,$(LIBS),-l$(LIBRARY)) -o $@

main_ind.o: main_ind.c

# Subset of eigenpairs by value
main_val.x: main_val.o
	$(LD) $(LDFLAGS) $< -L$(LIBPATH) \
        $(foreach LIBRARY,$(LIBS),-l$(LIBRARY)) -o $@

main_val.o: main_val.c

.PHONY: clean
clean:
	rm -f main_*.x core.* *__genmod.* *.o *~
