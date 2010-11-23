# Compiler for C and Fortran
CC = gcc
FC = gfortran

# Compiler flags
CFLAGS = -pthread -O3
FFLAGS = -O3 -funderscoring

# Archiver and flags used when building the archive
AR = /usr/bin/ar 
ARFLAGS = rcs

# Indicate if C99 feature of complex number are supported,
# which is ONLY NECESSARY FOR THE HERMITIAN WRAPPER ROUTINE.
# If true and routine needed, set to 1. To be safe it is set 
# to 0 by default
COMPLEX_SUPPORT = 0

# To build 'libmrrr.a' without adding necessary LAPACK routines 
# for the to the archive set value to 0; default value is 1
INCLAPACK = 1

# On some systems 'spinlocks' are not supported, therefore 
# here the flag to use 'mutexes' instead; default value is 1
SPINLOCK_SUPPORT = 1