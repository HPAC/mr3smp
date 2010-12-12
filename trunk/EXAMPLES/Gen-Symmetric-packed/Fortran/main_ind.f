      PROGRAM MAIN

      INTEGER N, NMAX, IL, IU, M, LDZ, SEED, SIZE
      PARAMETER (NMAX=1000, LDZ=NMAX,  
     $           SIZE=(NMAX*(NMAX+1))/2, SEED=1267893)

      DOUBLE PRECISION VL, VU, AP(SIZE), BP(SIZE), 
     $                 W(NMAX), Z(NMAX,NMAX)

      INTEGER ITYPE, I, J, K, IERR

*     external functions
      EXTERNAL DSPGEIG

*     Intialize symmetric matrices A and B (upper part stored)
      N = 100

      CALL SRAND(SEED)
      K = 1
      DO 100, J=1,N
         DO 200, I=1,J
            AP(K) = RAND()
            BP(K) = RAND()
            IF (I .EQ. J) THEN
               BP(K) = BP(K) + N
            ENDIF
            K = K + 1
 200     CONTINUE
 100  CONTINUE


*     Solve the symmetric eigenproblem with itype=1, that is 
*     A*x = lambda*B*x
*     The number of threads for the LAPACK routines are set by 
*     OMP_NUM_THREADS or GOTO_NUM_THREADS or MKL_NUM_THREADS ... 
*     depending on the BLAS used. For the tridiagonal stage with 
*     PMR_NUM_THREADS. 
      ITYPE = 1
      IL    = 1
      IU    = N/2 + 1
      CALL DSPGEIG(ITYPE, 'V', 'I', 'U', N, AP, BP,   
     $             VL, VU, IL, IU, M, W, Z, LDZ, IERR)
      IF (IERR .NE. 0) THEN
         WRITE(*,*) 'Routine has failed with error', IERR
      ENDIF


      WRITE(*,*) 'Sucessfully computed eigenpairs!'

      END
