      PROGRAM MAIN

      INTEGER N, NMAX, JMAX, IL, IU, M, LDA, LDZ, ZERO,
     $        SEED
      PARAMETER (NMAX=1000, JMAX=NMAX, LDA=NMAX, LDZ=NMAX, ZERO=0, 
     $           SEED=1)

      DOUBLE PRECISION VL, VU, A(NMAX,NMAX), W(NMAX), 
     $                 Z(NMAX,JMAX)

      INTEGER I, J, IERR

*     external functions
      EXTERNAL DSYEIG

*     Intialize symmetric matrix A of size N-by-N
      N = 300

      DO 100, I=1,N
         DO 200, J=1,I
            A(I,J) = RAND()
 200     CONTINUE
 100  CONTINUE

      DO 300, I=1,N
         DO 400, J=I+1,N
            A(I,J) = A(J,I)
 400     CONTINUE
 300  CONTINUE


*     Solve the symmetric eigenproblem
      CALL DSYEIG('V', 'A', 'L', N, A, LDA, VL, VU, IL, IU, 
     $            M, W, Z, LDZ, IERR)
      IF (IERR .NE. 0) THEN
         WRITE(*,*) 'Routine has failed with error', IERR
      ENDIF


      WRITE(*,*) 'Sucessfully computed eigenpairs!'

      END
