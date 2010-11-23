include ./make.inc

DIRS = ./INCLUDE ./SRC ./SRC/LAPACK ./SRC/LAPACK/SEQUENTIAL \
       ./SRC/BLAS

CFLAGS += -I./INCLUDE

ifeq ($(COMPLEX_SUPPORT),1)
 CFLAGS += -DCOMPLEX_SUPPORTED
endif

ifeq ($(SPINLOCK_SUPPORT),0)
 CFLAGS += -DNOSPINLOCKS
endif

# Source files
HEADERS := $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.h))
CSRCS   := $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.c))
COBJS = $(CSRCS:.c=.o)

ifeq ($(INCLAPACK),0)
 FSRCS   := $(foreach DIR,$(DIRS),$(wildcard $(DIR)/odscal.f))
 FOBJS = $(FSRCS:.f=.o)
else
 FSRCS   := $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.f))
 FOBJS = $(FSRCS:.f=.o)
endif

# Build target #
libmrrr.a: $(COBJS) $(FOBJS) $(HEADERS)
	   $(AR) $(ARFLAGS) ./LIB/libmrrr.a $(COBJS) $(FOBJS)

$(COBJS): $(HEADERS)
$(FOBJS):

.PHONY: clean
clean:
	rm -f *~ core.* *__genmod* \
        ./INSTALL/*~ \
        $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.o)) \
        $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*~)) \
        $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.mod.*)) \
        $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*__genmod*))