      LOGICAL FUNCTION ODNAN( DIN )
*
*  -- LAPACK auxiliary routine (version 3.2.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2010
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   DIN
*     ..
*
*  Purpose
*  =======
*
*  ODNAN returns .TRUE. if its argument is NaN, and .FALSE.
*  otherwise.  To be replaced by the Fortran 2003 intrinsic in the
*  future.
*
*  Arguments
*  =========
*
*  DIN     (input) DOUBLE PRECISION
*          Input to test for NaN.
*
*  =====================================================================
*
*  .. External Functions ..
      LOGICAL ODSNAN
      EXTERNAL ODSNAN
*  ..
*  .. Executable Statements ..
      ODNAN = ODSNAN(DIN,DIN)
      RETURN
      END
