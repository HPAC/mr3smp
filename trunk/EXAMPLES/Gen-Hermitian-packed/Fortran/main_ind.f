      PROGRAM MAIN

      INTEGER N, NMAX, IL, IU, M, LDZ, ZERO, SIZE, SEED
      PARAMETER (NMAX=1000, LDZ=NMAX, ZERO=0, 
     $           SIZE=(NMAX*(NMAX+1))/2, SEED=976733)

      DOUBLE PRECISION VL, VU, W(NMAX)      
      COMPLEX*16       AP(SIZE), BP(SIZE), Z(NMAX,NMAX)

      INTEGER ITYPE, I, J, K, IERR

*     external functions
      EXTERNAL ZHPGEIG

*     Intialize Hermitian matrix A of size N-by-N
      N = 100

      CALL SRAND(SEED)
      K = 1
      DO 100, J=1,N
         DO 200, I=1,J
            IF (I .EQ. J) THEN
               AP(K) = COMPLEX(RAND(),ZERO)
               BP(K) = COMPLEX(RAND(),ZERO) + N
            ELSE
               AP(K) = COMPLEX(RAND(),RAND())
               BP(K) = COMPLEX(RAND(),RAND())
            ENDIF
            K = K + 1
 200     CONTINUE
 100  CONTINUE


*     Solve the eigenproblem of for A*x = lambda*B*x
*     The number of threads for the LAPACK routines are set by 
*     OMP_NUM_THREADS or GOTO_NUM_THREADS or MKL_NUM_THREADS ... 
*     depending on the BLAS used. For the tridiagonal stage with 
*     PMR_NUM_THREADS. 
      ITYPE = 1
      IL    = 2
      IU    = 3
      CALL ZHPGEIG(ITYPE, 'V', 'I', 'U', N, AP, BP,   
     $             VL, VU, IL, IU, M, W, Z, LDZ, IERR)
      IF (IERR .NE. 0) THEN
         WRITE(*,*) 'Routine has failed with error', IERR
      ENDIF


      WRITE(*,*) 'Sucessfully computed eigenpairs!'

      END
