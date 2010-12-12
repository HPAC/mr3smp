      PROGRAM MAIN

      INTEGER N, NMAX, JMAX, IL, IU, M, LDZ, ZERO,
     $        SIZE, SEED
      PARAMETER (NMAX=1000, JMAX=NMAX, LDZ=NMAX, ZERO=0, 
     $           SIZE=(NMAX*(NMAX+1))/2, SEED=7873)

      DOUBLE PRECISION VL, VU, W(NMAX)      
      COMPLEX*16       AP(SIZE), Z(NMAX,JMAX)

      INTEGER I, J, K, IERR

*     external functions
      EXTERNAL ZHPEIG

*     Intialize Hermitian matrix A of size N-by-N, where the 
*     storage of the upper part is assumed
      N = 100

      CALL SRAND(SEED)
      K = 1
      DO 100, J=1,N
         DO 200, I=1,J
            IF (I .EQ. J) THEN
               AP(K) = COMPLEX(RAND(),ZERO)
            ELSE
               AP(K) = COMPLEX(RAND(),RAND())
            ENDIF
            K =  K + 1
 200     CONTINUE
 100  CONTINUE


*     Solve the eigenproblem
*     The number of threads for the LAPACK routines are set by 
*     OMP_NUM_THREADS or GOTO_NUM_THREADS or MKL_NUM_THREADS ... 
*     depending on the BLAS used. For the tridiagonal stage with 
*     PMR_NUM_THREADS. 
      IL = 1
      IU = N/2 + 1
      CALL ZHPEIG('V', 'I', 'U', N, AP, VL, VU, IL, IU, 
     $            M, W, Z, LDZ, IERR)
      IF (IERR .NE. 0) THEN
         WRITE(*,*) 'Routine has failed with error', IERR
      ENDIF


      WRITE(*,*) 'Sucessfully computed eigenpairs!'

      END
