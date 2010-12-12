      PROGRAM MAIN

      INTEGER N, NMAX, JMAX, IL, IU, M, LDA, LDZ, SEED
      PARAMETER (NMAX=1000, JMAX=NMAX, LDA=NMAX, LDZ=NMAX,  
     $           SEED=13)

      DOUBLE PRECISION VL, VU, A(NMAX,NMAX), W(NMAX), 
     $                 Z(NMAX,JMAX)

      INTEGER I, J, IERR

*     external functions
      EXTERNAL DSYEIG

*     Intialize symmetric matrix A of size N-by-N
      N = 100

      CALL SRAND(SEED)
      DO 100, J=1,N
         DO 200, I=1,J
            A(I,J) = RAND()
 200     CONTINUE
 100  CONTINUE

      DO 300, J=1,N
         DO 400, I=J+1,N
            A(I,J) = A(J,I)
 400     CONTINUE
 300  CONTINUE


*     Solve the symmetric eigenproblem
*     The number of threads for the LAPACK routines are set by 
*     OMP_NUM_THREADS or GOTO_NUM_THREADS or MKL_NUM_THREADS ... 
*     depending on the BLAS used. For the tridiagonal stage with 
*     PMR_NUM_THREADS. 
      IL = 1
      IU = N/2 + 1
      CALL DSYEIG('V', 'I', 'L', N, A, LDA, VL, VU, IL, IU, 
     $            M, W, Z, LDZ, IERR)
      IF (IERR .NE. 0) THEN
         WRITE(*,*) 'Routine has failed with error', IERR
      ENDIF


      WRITE(*,*) 'Sucessfully computed eigenpairs!'

      END
