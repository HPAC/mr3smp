      PROGRAM MAIN

      INTEGER N, NMAX, IL, IU, M, LDA, LDB, ZERO,
     $        SEED
      PARAMETER (NMAX=1000, LDA=NMAX, LDB=NMAX, ZERO=0, 
     $           SEED=7673263)

      DOUBLE PRECISION VL, VU, W(NMAX)      
      COMPLEX*16 A(NMAX,NMAX), B(NMAX,NMAX)

      INTEGER ITYPE, I, J, IERR

*     external functions
      EXTERNAL ZHEGEIG

*     Intialize Hermitian matrix A of size N-by-N
      N = 100

      CALL SRAND(SEED)
      DO 100, J=1,N
         DO 200, I=1,J
            IF (I .EQ. J) THEN
               A(I,J) = COMPLEX(RAND(),ZERO)
               B(I,J) = COMPLEX(RAND(),ZERO) + N
            ELSE
               A(I,J) = COMPLEX(RAND(),RAND())
               B(I,J) = COMPLEX(RAND(),RAND())
            ENDIF
 200     CONTINUE
 100  CONTINUE

      DO 300, J=1,N
         DO 400, I=J+1,N
            A(I,J) = CONJG(A(J,I))
            B(I,J) = CONJG(B(J,I))
 400     CONTINUE
 300  CONTINUE



*     Solve the eigenproblem of for A*x = lambda*B*x
*     The number of threads for the LAPACK routines are set by 
*     OMP_NUM_THREADS or GOTO_NUM_THREADS or MKL_NUM_THREADS ... 
*     depending on the BLAS used. For the tridiagonal stage with 
*     PMR_NUM_THREADS. 
      ITYPE=1
      CALL ZHEGEIG(ITYPE, 'V', 'A', 'L', N, A, LDA, B, LDB,  
     $             VL, VU, IL, IU, M, W, IERR)
      IF (IERR .NE. 0) THEN
         WRITE(*,*) 'Routine has failed with error', IERR
      ENDIF



      WRITE(*,*) 'Sucessfully computed eigenpairs!'

      END
