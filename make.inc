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

# Compile routines without adding necessary LAPACK routines 
# to 'libmrrr.a'; default = 0
NOLAPACK = 0