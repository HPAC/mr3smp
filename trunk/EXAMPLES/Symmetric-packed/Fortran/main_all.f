      PROGRAM MAIN

      INTEGER N, NMAX, JMAX, IL, IU, M, LDA, LDZ, ZERO, SEED,
     $        SIZE
      PARAMETER (NMAX=1000, JMAX=NMAX, LDA=NMAX, LDZ=NMAX, ZERO=0, 
     $           SIZE=(NMAX*(NMAX+1))/2, SEED=13)

      DOUBLE PRECISION VL, VU, AP(SIZE), W(NMAX), Z(NMAX,JMAX)

      INTEGER I, J, K, IERR

*     external functions
      EXTERNAL DSPEIG

*     Intialize symmetric matrix A of size N-by-N
      N = 100

      CALL SRAND(SEED)
      K = 1
      DO 100, J=1,N
         DO 200, I=1,J
            AP(K) = RAND()
            K     = K + 1
 200     CONTINUE
 100  CONTINUE


*     Solve the eigenproblem
*     The number of threads for the LAPACK routines are set by 
*     OMP_NUM_THREADS or GOTO_NUM_THREADS or MKL_NUM_THREADS ... 
*     depending on the BLAS used. For the tridiagonal stage with 
*     PMR_NUM_THREADS. 
      CALL DSPEIG('V', 'A', 'U', N, AP, VL, VU, IL, IU, 
     $            M, W, Z, LDZ, IERR)
      IF (IERR .NE. 0) THEN
         WRITE(*,*) 'Routine has failed with error', IERR
      ENDIF


      WRITE(*,*) 'Sucessfully computed eigenpairs!'

      END
