      PROGRAM MAINA

      INTEGER N, NMAX, JMAX, IL, IU, TRYRAC, M, LDZ, ZERO
      PARAMETER (NMAX=1000, JMAX=NMAX, LDZ=NMAX, ZERO=0)

      DOUBLE PRECISION VL, VU, D(NMAX), E(NMAX), W(NMAX), 
     $                 Z(NMAX,JMAX), ZSUPP(2*NMAX)

      INTEGER I, J, UFILE, IERR
      PARAMETER (UFILE=1)

*     external functions
      EXTERNAL MRRR

*     Read in data
      OPEN(UFILE, FILE='Wilkinson21.data', STATUS='OLD', IOSTAT=IERR,
     $     FORM='FORMATTED')

      READ(UFILE, 10) N
      DO 100, I=1,N
         READ(UFILE, 20, END=200) D(I), E(I)
 100  CONTINUE
 200  CONTINUE

      CLOSE(UFILE)

*     Call MRRR ( if underscores not added automatically MRRR_() )
      TRYRAC = 1

      CALL MRRR('V', 'A', N, D, E, VL, VU, IL, IU, TRYRAC, 
     $          M, W, Z, LDZ, ZSUPP, IERR)
      IF (IERR .NE. 0) THEN
         WRITE(*,*) 'Routine has failed with error', IERR
      ENDIF

*     Write out results; this is just a hack and will sometimes fail
*     since the output different processes may interleave
      OPEN(UFILE, FILE='result_all.m', STATUS='UNKNOWN')
      CLOSE(UFILE)

      OPEN(UFILE, FILE='result_all.m', STATUS='OLD', 
     $     ACCESS='APPEND')
      CALL PRTVEC(UFILE, 'W', W, M, ZERO)
      
      OPEN(UFILE, FILE='result_all.m', STATUS='OLD', 
     $     ACCESS='APPEND')
      DO 500, J=1,M
         CALL PRTCOL(UFILE, 'Z', Z(1,J), N, J, ZERO)
 500  CONTINUE
      CLOSE(UFILE)

 10   FORMAT(1X, I4)
 20   FORMAT(1X, E23.17, 1X, E23.17E2)

      END




      SUBROUTINE PRTVEC(UFILE, NAME, V, N, OFFSET)

      CHARACTER NAME
      INTEGER I, N, OFFSET, UFILE
      DOUBLE PRECISION V(N)

      DO 100, I=1,N
         WRITE(UFILE,10) NAME,'(',I+OFFSET,')=',V(I),';'
 100  CONTINUE

 10   FORMAT(1X, A1, A1, I4, A2, 1X, E25.17E3, A1)
      
      END




      SUBROUTINE PRTCOL(UFILE, NAME, V, N, COL, OFFSET)

      CHARACTER NAME
      INTEGER I, N, COL, OFFSET, UFILE
      DOUBLE PRECISION V(N)

      DO 100, I=1,N
         WRITE(UFILE,10) NAME,'(',I,',',COL+OFFSET,')=',V(I),';'
 100  CONTINUE

 10   FORMAT(1X, A1, A1, I4, A1, I4, A2, 1X, E25.17E3, A1)
      
      END
