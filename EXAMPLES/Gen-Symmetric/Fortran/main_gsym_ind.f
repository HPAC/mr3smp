      PROGRAM MAIN

      INTEGER N, NMAX, IL, IU, M, LDA, LDB, ZERO,
     $        SEED
      PARAMETER (NMAX=1000, LDA=NMAX, LDB=NMAX, ZERO=0, 
     $           SEED=1267893)

      DOUBLE PRECISION VL, VU, A(NMAX,NMAX), B(NMAX,NMAX), 
     $                 W(NMAX)

      INTEGER ITYPE, I, J, IERR

*     external functions
      EXTERNAL DSYGEIG

*     Intialize symmetric matrices A and B of size N-by-N full
*     such that either only the upper or lower triangular part 
*     can later be used
      N = 100

      CALL SRAND(SEED)
      DO 100, J=1,N
         DO 200, I=1,J
            A(I,J) = RAND()
            B(I,J) = RAND()
            IF (I .EQ. J) THEN
               B(I,J) = B(I,J) + N
            ENDIF
 200     CONTINUE
 100  CONTINUE

      DO 300, J=1,N
         DO 400, I=J+1,N
            A(I,J) = A(J,I)
            B(I,J) = B(J,I)
 400     CONTINUE
 300  CONTINUE



*     Solve the symmetric eigenproblem with itype=1, that is 
*     A*x = lambda*B*x
*     Array A contains the eigenvectors on exit
*     The number of threads for the LAPACK routines are set by 
*     OMP_NUM_THREADS or GOTO_NUM_THREADS or MKL_NUM_THREADS ... 
*     depending on the BLAS used. For the tridiagonal stage with 
*     PMR_NUM_THREADS. 
      ITYPE = 1
      IL    = 1
      IU    = N/2 + 1
      CALL DSYGEIG(ITYPE, 'V', 'I', 'L', N, A, LDA, B, LDB,  
     $             VL, VU, IL, IU, M, W, IERR)
      IF (IERR .NE. 0) THEN
         WRITE(*,*) 'Routine has failed with error', IERR
      ENDIF



      WRITE(*,*) 'Sucessfully computed eigenpairs!'

      END
