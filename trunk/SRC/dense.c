/* Routines for the dense symmetric and hermitian eigenproblem 
 * using the multithreaded MRRR for the tridiagonal stage.
 *
 * See README and header file 'mrrr.h' for more information.
 *
 * Copyright (c) 2010, RWTH Aachen University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or 
 * without modification, are permitted provided that the following
 * conditions are met:
 *   * Redistributions of source code must retain the above 
 *     copyright notice, this list of conditions and the following
 *     disclaimer.
 *   * Redistributions in binary form must reproduce the above 
 *     copyright notice, this list of conditions and the following 
 *     disclaimer in the documentation and/or other materials 
 *     provided with the distribution.
 *   * Neither the name of the RWTH Aachen University nor the
 *     names of its contributors may be used to endorse or promote 
 *     products derived from this software without specific prior 
 *     written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL RWTH 
 * AACHEN UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF 
 * SUCH DAMAGE.
 *
 * Coded by Matthias Petschow (petschow@aices.rwth-aachen.de),
 * September 2010, Version 1.2
 *
 * This code was the result of a collaboration between 
 * Matthias Petschow and Paolo Bientinesi. When you use this 
 * code, kindly reference the paper:
 *
 * "MR3-SMP: A Symmetric Tridiagonal Eigensolver for Multi-Core 
 * Architectures" by Matthias Petschow and Paolo Bientinesi, 
 * RWTH Aachen, Technical Report AICES-2010/10-2 (submitted to 
 * Parallel Computing).
 *
 */


#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <assert.h>
#ifdef COMPLEX_SUPPORTED
#include <complex.h>
#endif
#include "global.h"
#include "mrrr.h"  

/* LAPACK routines for functions of the symmetric and the Hermitian 
 * eigenproblem */ 
extern void dsytrd_(char*, int*, double*, int*, double*, double*, 
		    double*, double*, int*, int*);
extern void dormtr_(char*, char*, char*, int*, int*, double*, int*, 
		    double*, double*, int*, double*, int*, int*);
extern double dlansy_(char*, char*, int*, double*, int*, double*);
#ifdef COMPLEX_SUPPORTED
extern void zhetrd_(char*, int*, double complex*, int*, double*, 
		    double*, double complex*, double complex*, int*, 
		    int*);
extern void zunmtr_(char*, char*, char*, int*, int*, double complex*, 
		    int*, double complex*, double complex*, int*, 
		    double complex*, int*, int*);
extern double zlanhe_(char*, char*, int*, double complex*, int*, 
		      double complex*);
extern void zdscal_(int*, double*, double complex*, int*);
#endif

/* LAPACK routines for functions of the generalized symmetric-definite 
 * and the Hermitian-definite eigenproblem */ 
extern void dpotrf_(char*, int*, double*, int*, int*);
extern void dsygst_(int*, char*, int*, double*, int*, double*, 
		    int*, int*);
extern void dtrsm_(char*, char*, char*, char*, int*, int*, double*, 
		   double*, int*, double*, int*);
extern void dtrmm_(char*, char*, char*, char*, int*, int*, double*, 
		   double*, int*, double*, int*);
#ifdef COMPLEX_SUPPORTED
extern void zpotrf_(char*, int*, double complex*, int*, int*);
extern void zhegst_(int*, char*, int*, double complex*, int*, 
		    double complex*, int*, int*);
extern void ztrsm_(char*, char*, char*, char*, int*, int*, 
		   double complex*, double complex*, int*, 
		   double complex*, int*);
extern void ztrmm_(char*, char*, char*, char*, int*, int*, 
		   double complex*, double complex*, int*,
		   double complex*, int*);
#endif

/* LAPACK routines for functions of the symmetric and the Hermitian 
 * eigenproblem in packed storage */ 
extern void dsptrd_(char*,int*,double*,double*,double*,double*,int*);
extern void dopmtr_(char*,char*,char*,int*,int*,double*,double*,
		    double*,int*,double*,int*);
extern double dlansp_(char*,char*,int*,double*,double*);
#ifdef COMPLEX_SUPPORTED
extern void zhptrd_(char*,int*,double complex*,double*,double*,
		    double complex*,int*);
extern void zupmtr_(char*,char*,char*,int*,int*,double complex*,
		    double complex*,double complex*,int*,
		    double complex*,int*);
extern double zlanhp_(char*,char*,int*,double complex*,double*);
#endif

/* LAPACK routines for functions of the generalized symmetric-definite 
 * and the Hermitian-definite eigenproblem in packed storage */
extern void dpptrf_(char*,int*,double*,int*);
extern void dspgst_(int*,char*,int*,double*,double*,int*);
extern void dtpsv_(char*,char*,char*,int*,double*,double*,int*);
extern void dtpmv_(char*,char*,char*,int*,double*,double*,int*); 
#ifdef COMPLEX_SUPPORTED
extern void zpptrf_(char*,int*,double complex*,int*);
extern void zhpgst_(int*,char*,int*,double complex*,
		    double complex*,int*);
extern void ztpsv_(char*,char*,char*,int*,double complex*,
		   double complex*,int*);
extern void ztpmv_(char*,char*,char*,int*,double complex*,
		   double complex*,int*);
#endif


/* Other prototypes */
static double dscale_matrix(char*, char*, int*, double*, int*, 
			    double*, double*, double*);
static double dscale_matrix_packed(char*, char*, int*, double*, 
				   double*, double*, double*);
#ifdef COMPLEX_SUPPORTED
static double zscale_matrix(char*, char*, int*, double complex*, int*,
			    double*, double*, double complex*);
static double zscale_matrix_packed(char*, char*, int*, 
				   double complex*, double*, 
				   double*, double complex*);
#endif



/* Routine for the dense symmetric eigenproblem */
int dsyeig(char *jobz, char *range, char *uplo, int *np, double *A, 
	   int *ldap, double *vlp, double *vup, int *ilp, int *iup, 
	   int *mp, double *W, double *Z, int *ldzp)
{
  int    n      = *np;
  int    szwrk  = 512*n;
  int    tryRAC = 1;

  double *D, *E, *TAU;
  double *work;
  int    *Zsupp;
  int    info, ione=1;

  double scale = 1.0;
  double invscale;

  bool   onlyW  = (jobz[0]  == 'N' || jobz[0]  == 'n');
  bool   wantZ  = (jobz[0]  == 'V' || jobz[0]  == 'v');
  bool   cntval = (jobz[0]  == 'C' || jobz[0]  == 'c');
  bool   alleig = (range[0] == 'A' || range[0] == 'a');
  bool   valeig = (range[0] == 'V' || range[0] == 'v');
  bool   indeig = (range[0] == 'I' || range[0] == 'i');

  if( !(onlyW  || wantZ  || cntval) ) return(1);
  if( !(alleig || valeig || indeig) ) return(1);
  if(n <= 0) return(1);
  if (valeig) {
    if(*vup<=*vlp) return(1);
  } else if (indeig) {
    if (*ilp<1 || *ilp>n || *iup<*ilp || *iup>n) return(1);
  }

  D = (double *) malloc( n*sizeof(double) );
  assert(D != NULL);

  E = (double *) malloc( n*sizeof(double) );
  assert(E != NULL);

  TAU = (double *) malloc( n*sizeof(double) );
  assert(TAU != NULL);

  work = (double *) malloc( szwrk*sizeof(double) );
  assert(work != NULL);

  Zsupp = (int *) malloc( 2*n*sizeof(int) );
  assert(Zsupp != NULL);

  /* Scale matrix if necessary */
  scale = dscale_matrix(range, uplo, np, A, ldap, vlp, vup, work);

  /* Reduction to tridiagonal */
  dsytrd_(uplo, np, A, ldap, D, E, TAU, work, &szwrk, &info);
  assert(info == 0);

  /* Use MRRR to compute eigenvalues and -vectors */
  info = mrrr(jobz, range, np, D, E, vlp, vup, ilp, iup, 
	      &tryRAC, mp, W, Z, ldzp, Zsupp);
  assert(info == 0);
  
  /* Backtransformation Z = Q*Z */
  dormtr_("L", uplo, "N", np, mp, A, ldap, TAU, Z, ldzp, work,
	  &szwrk, &info);
  assert(info == 0);

  /* Scaling of eigenvalues if necessary */
  if (scale != 1.0) { /* FP cmp okay */
    invscale = 1.0/scale;
    *vlp *= invscale;
    *vup *= invscale;
    dscal_(mp, &invscale, W, &ione);
  }

  free(D);
  free(E);
  free(TAU);
  free(work);
  free(Zsupp);

  return(0);
}




/* Routine for the dense generalized symmetric-definite eigenproblem */
int dsygeig(int *itype, char *jobz, char *range, char *uplo, int *np, 
	    double *A, int *ldap, double *B, int *ldbp, double *vlp, 
	    double *vup, int *ilp, int *iup, int *mp, double *W)
{
  int    info, itmp, j;
  bool   onlyW  = (jobz[0]  == 'N' || jobz[0]  == 'n');
  bool   wantZ  = (jobz[0]  == 'V' || jobz[0]  == 'v');
  bool   cntval = (jobz[0]  == 'C' || jobz[0]  == 'c');
  bool   upper  = (uplo[0]  == 'U' || uplo[0]  == 'u');
  bool   lower  = (uplo[0]  == 'L' || uplo[0]  == 'l');
  bool   alleig = (range[0] == 'A' || range[0] == 'a');
  bool   valeig = (range[0] == 'V' || range[0] == 'v');
  bool   indeig = (range[0] == 'I' || range[0] == 'i');
  double one    = 1.0;
  int    n      = *np;
  int    lda    = *ldap;
  char   *trans;
  double *Z;

  /* Check input */
  if (n <= 1) return(1);
  if (*itype < 1 || *itype > 3) return(1);
  if (!lower  && !upper) return(1);
  if (!alleig && !valeig && !indeig) return(1);
  if (!onlyW  && !wantZ  && !cntval) return(1);

  if (indeig) itmp = *iup-*ilp+1;
  else        itmp = n;
  Z = (double *) malloc((size_t) itmp*n*sizeof(double));
  assert(Z != NULL);

  /* Form the Cholesky factor of B */
  dpotrf_(uplo, np, B, ldbp, &info);
  assert(info == 0);

  /* Convert problem to standard eigenvalue problem */
  dsygst_(itype, uplo, np, A, ldap, B, ldbp, &info);
  assert(info == 0);

  /* Solve standard eigenvalue problem using MRRR */
  info = dsyeig(jobz, range, uplo, np, A, ldap, vlp, vup, 
		ilp, iup, mp, W, Z, np);
  assert(info == 0);

  /* Backtransform eigenvectors */
  if (wantZ) {

    if (*itype ==  1 || *itype == 2) {
      /* A*x = lambda*B*x or A*B*x = lambda*x requires 
       * x = inv(L)'*y or x = inv(U)*y */

      if (upper) trans = "N";
      else       trans = "T";

      dtrsm_("Left", uplo, trans, "Non-unit", np, mp, &one, 
	     B, ldbp, Z, np);
    } else if (*itype == 3) {
      /* B*A*x = lambda*x requires x = L*y or U'*y */
      
      if (upper) trans = "T";
      else       trans = "N";

      dtrmm_("Left", uplo, trans, "Non-unit", np, mp, &one, 
	     B, ldbp, Z, np);
    } else {
      return(1);
    }

  } /* Backtransform eigenvectors */ 

  /* Copy eigenvectors into A */
  for (j=0; j<*mp; j++)
    memcpy(&A[j*lda], &Z[j*n], n*sizeof(double));

  free(Z);

  return(0);
}




/* Routine for the dense symmetric eigenproblem in packed storage */
int dspeig(char *jobz, char *range, char *uplo, int *np, double *AP, 
	   double *vlp, double *vup, int *ilp, int *iup, 
	   int *mp, double *W, double *Z, int *ldzp)
{
  int    n      = *np;
  int    tryRAC = 1;

  double *D, *E, *TAU;
  double *work;
  int    *Zsupp;
  int    info, ione=1;

  double scale = 1.0;
  double invscale;

  bool   onlyW  = (jobz[0]  == 'N' || jobz[0]  == 'n');
  bool   wantZ  = (jobz[0]  == 'V' || jobz[0]  == 'v');
  bool   cntval = (jobz[0]  == 'C' || jobz[0]  == 'c');
  bool   alleig = (range[0] == 'A' || range[0] == 'a');
  bool   valeig = (range[0] == 'V' || range[0] == 'v');
  bool   indeig = (range[0] == 'I' || range[0] == 'i');

  if( !(onlyW  || wantZ  || cntval) ) return(1);
  if( !(alleig || valeig || indeig) ) return(1);
  if(n <= 0) return(1);
  if (valeig) {
    if(*vup<=*vlp) return(1);
  } else if (indeig) {
    if (*ilp<1 || *ilp>n || *iup<*ilp || *iup>n) return(1);
  }
  
  D = (double *) malloc( n*sizeof(double) );
  assert(D != NULL);
  
  E = (double *) malloc( n*sizeof(double) );
  assert(E != NULL);
  
  TAU = (double *) malloc( n*sizeof(double) );
  assert(TAU != NULL);
  
  work = (double *) malloc( n*sizeof(double) );
  assert(work != NULL);
  
  Zsupp = (int *) malloc( 2*n*sizeof(int) );
  assert(Zsupp != NULL);
  
  /* Scale matrix if necessary */
  scale = dscale_matrix_packed(range, uplo, np, AP, vlp, vup, work);

  /* Reduction to tridiagonal */
  dsptrd_(uplo, np, AP, D, E, TAU, &info);
  assert(info == 0);
  
  /* Use MRRR to compute eigenvalues and -vectors */
  info = mrrr(jobz, range, np, D, E, vlp, vup, ilp, iup, 
  	      &tryRAC, mp, W, Z, ldzp, Zsupp);
  assert(info == 0);
  
  /* Backtransformation Z = Q*Z */
  dopmtr_("L", uplo, "N", np, mp, AP, TAU, Z, ldzp, work, &info);
  assert(info == 0);

  /* Scaling of eigenvalues if necessary */
  if (scale != 1.0) { /* FP cmp okay */
    invscale = 1.0/scale;
    *vlp *= invscale;
    *vup *= invscale;
    dscal_(mp, &invscale, W, &ione);
  }
  
  free(D);
  free(E);
  free(TAU);
  free(work);
  free(Zsupp);
  
  return(0);
}




/* Routine for the dense generalized symmetric-definite eigenproblem 
 * using packed storage */
int dspgeig(int *itype, char *jobz, char *range, char *uplo, int *np, 
	    double *AP, double *BP, double *vlp, double *vup, 
	    int *ilp, int *iup, int *mp, double *W, double *Z, 
	    int *ldzp)
{
  int    n      = *np;
  int    ldz    = *ldzp;
  int    info, j, ione=1;
  bool   onlyW  = (jobz[0]  == 'N' || jobz[0]  == 'n');
  bool   wantZ  = (jobz[0]  == 'V' || jobz[0]  == 'v');
  bool   cntval = (jobz[0]  == 'C' || jobz[0]  == 'c');
  bool   upper  = (uplo[0]  == 'U' || uplo[0]  == 'u');
  bool   lower  = (uplo[0]  == 'L' || uplo[0]  == 'l');
  bool   alleig = (range[0] == 'A' || range[0] == 'a');
  bool   valeig = (range[0] == 'V' || range[0] == 'v');
  bool   indeig = (range[0] == 'I' || range[0] == 'i');
  char   *trans;

  /* Check input */
  if (n < 1) return(1);
  if (*itype < 1 || *itype > 3) return(1);
  if (!lower  && !upper) return(1);
  if (!alleig && !valeig && !indeig) return(1);
  if (!onlyW  && !wantZ  && !cntval) return(1);

  /* Form the Cholesky factor of B */
  dpptrf_(uplo, np, BP, &info);
  assert(info == 0);

  /* Convert problem to standard eigenvalue problem */
  dspgst_(itype, uplo, np, AP, BP, &info);
  assert(info == 0);

  /* Solve standard eigenvalue problem using MRRR */
  info = dspeig(jobz, range, uplo, np, AP, vlp, vup, 
  		ilp, iup, mp, W, Z, ldzp);
  assert(info == 0);

  /* Backtransform eigenvectors */
  if (wantZ) {

    if (*itype ==  1 || *itype == 2) {
      /* A*x = lambda*B*x or A*B*x = lambda*x requires 
       * x = inv(L)'*y or x = inv(U)*y */

      if (upper) trans = "N";
      else       trans = "T";

      for (j=0; j<*mp; j++)
	dtpsv_(uplo, trans, "Non-unit", np, BP, &Z[j*ldz], &ione);

    } else if (*itype == 3) {
      /* B*A*x = lambda*x requires x = L*y or U'*y */
  
      if (upper) trans = "T";
      else       trans = "N";

      for (j=0; j<*mp; j++)
	dtpmv_(uplo, trans, "Non-unit", np, BP, &Z[j*ldz], &ione);

    } else {
      return(1);
    }
  
  } /* Backtransform eigenvectors */ 

  return(0);
}

    


/* Scale dense matrix to allowable range */
static 
double dscale_matrix(char *range, char *uplo, int *np, double *A, 
		     int *ldap, double *vlp, double *vup, double *work)
{
  double sigma = 1.0;
  double smlnum, bignum, rmin, rmax;
  double norm;
  bool   scaled = false;
  bool   valeig = (range[0] == 'V' || range[0] == 'v');
  bool   lower  = (uplo[0] == 'L' || uplo[0] == 'l');
  int    n      = *np;
  int    lda    = *ldap;
  int    i, ione=1, itmp;

  smlnum = DBL_MIN / DBL_EPSILON;
  bignum = 1.0 / smlnum;
  rmin   = sqrt(smlnum);
  rmax   = fmin(sqrt(bignum), 1.0 / sqrt(sqrt(DBL_MIN)));

  norm = dlansy_("M", uplo, np, A, ldap, work);
  if (norm > 0.0 && norm < rmin) {
    scaled = true;
    sigma  = rmin / norm;
  } else if (norm > rmax) {
    scaled = true;
    sigma  = rmax / norm;
  }
  if (scaled) {
    if (lower) {
      for (i=0; i<n; i++) {
	itmp = n - i;
	dscal_(&itmp, &sigma, &A[i + i*lda], &ione);
      }
    } else {
      for (i=0; i<n; i++) {
	itmp = i + 1;
	dscal_(&itmp, &sigma, &A[i*lda], &ione);
      }
    }
    if (valeig) {
      *vlp *= sigma;
      *vup *= sigma;
    }
  }
  
  return(sigma);
}
  


  
/* Scale matrix in packed storage to allowable range */ 
static 
double dscale_matrix_packed(char *range, char *uplo, int *np, double *AP, 
			    double *vlp, double *vup, double *work)
{
  double sigma = 1.0;
  double smlnum, bignum, rmin, rmax;
  double norm;
  bool   scaled = false;
  bool   valeig = (range[0] == 'V' || range[0] == 'v');
  int    n      = *np;
  int    ione=1, itmp;

  smlnum = DBL_MIN / DBL_EPSILON;
  bignum = 1.0 / smlnum;
  rmin   = sqrt(smlnum);
  rmax   = fmin(sqrt(bignum), 1.0 / sqrt(sqrt(DBL_MIN)));

  norm = dlansp_("M", uplo, np, AP, work);
  if (norm > 0.0 && norm < rmin) {
    scaled = true;
    sigma  = rmin / norm;
  } else if (norm > rmax) {
    scaled = true;
    sigma  = rmax / norm;
  }
  if (scaled) {
    itmp = (n*(n+1))/2;
    dscal_(&itmp, &sigma, AP, &ione);
    if (valeig) {
      *vlp *= sigma;
      *vup *= sigma;
    }
  }
  
  return(sigma);
}



#ifdef COMPLEX_SUPPORTED
int zheeig(char *jobz, char *range, char *uplo, int *np, 
	   double complex *A, int *ldap, double *vlp, double *vup, 
	   int *ilp, int *iup, int *mp, double *W, double complex *Z, 
	   int *ldzp)
{
  int n       = *np;
  long int nn = n;
  int szwrk   = 512*n;
  int tryRAC  = 1;

  double         *D, *E, *Ztmp;
  double complex *TAU;
  double complex *work;
  int            *Zsupp;
  long int       tmp, mm;
  long int       i, j;
  int            m, info;
  double         scale = 1.0;
  double         invscale;
  int            ione=1;
  
  bool   onlyW  = (jobz[0]  == 'N' || jobz[0]  == 'n');
  bool   wantZ  = (jobz[0]  == 'V' || jobz[0]  == 'v');
  bool   cntval = (jobz[0]  == 'C' || jobz[0]  == 'c');
  bool   alleig = (range[0] == 'A' || range[0] == 'a');
  bool   valeig = (range[0] == 'V' || range[0] == 'v');
  bool   indeig = (range[0] == 'I' || range[0] == 'i');

  if( !(onlyW  || wantZ  || cntval) ) return(1);
  if( !(alleig || valeig || indeig) ) return(1);
  if(n <= 0) return(1);
  if (valeig) {
    if(*vup<=*vlp) return(1);
  } else if (indeig) {
    if (*ilp<1 || *ilp>n || *iup<*ilp || *iup>n) return(1);
  }

  D = (double *) malloc(n*sizeof(double));
  assert(D != NULL);

  E = (double *) malloc(n*sizeof(double));
  assert(E != NULL);

  TAU = (double complex *) malloc((size_t) n*sizeof(double complex));
  assert(TAU != NULL);

  work = (double complex *) malloc( szwrk*sizeof(double complex) );
  assert(work != NULL);

  Zsupp = (int *) malloc(2*n*sizeof(int));
  assert(Zsupp != NULL);

  /* Scale matrix if necessary */
  scale = zscale_matrix(range, uplo, np, A, ldap, vlp, vup, work);

  /* Reduction to tridiagonal */
  zhetrd_(uplo, np, A, ldap, D, E, TAU, work, &szwrk, &info);
  assert(info == 0);

  /* Use MRRR to compute eigenvalues and -vectors, where part of 
   * Z (starting at Ztmp) is used to store the real eigenvectors 
   * of the tridiagonal temporarily */
  if (alleig)
    mm    = n;
  else if (indeig)
    mm    = (*iup)-(*ilp)+1;
  else {
    info = mrrr("Count", range, np, D, E, vlp, vup, ilp, iup, 
		&tryRAC, &m, W, NULL, ldzp, Zsupp);
    mm = m;
  }
  tmp  = (nn*mm)/2 + ( ((nn*mm) % 2) > 0 ); /* ceil(n*m/2) */
  Ztmp = (double *) &Z[(*ldzp) * mm - tmp];

  info = mrrr(jobz, range, np, D, E, vlp, vup, ilp, iup, 
	      &tryRAC, mp, W, Ztmp, np, Zsupp);
  assert(info == 0);

  /* Copy intermediate real eigenvectors to complex Z */ 
  for (j=0; j<(*mp); j++)
    for (i=0; i<(*np); i++)
      Z[j * (*ldzp) + i] = Ztmp[j * (*np) + i] + 0.0 * I;
  
  /* Backtransformation Z = U*Z */
  zunmtr_("L", uplo, "N", np, mp, A, ldap, TAU, Z, ldzp, 
	  work, &szwrk, &info);
  assert(info == 0);

  /* Scaling of eigenvalues if necessary */
  if (scale != 1.0) { /* FP cmp okay */
    invscale = 1.0/scale;
    *vlp *= invscale;
    *vup *= invscale;
    dscal_(mp, &invscale, W, &ione);
  }

  free(D);
  free(E);
  free(TAU);
  free(work);
  free(Zsupp);

  return(0);
}



/* Routine for the dense generalized Hermitian-definite eigenproblem */
int zhegeig(int *itype, char *jobz, char *range, char *uplo, int *np, 
	    double complex *A, int *ldap, double complex *B, int *ldbp, 
	    double *vlp, double *vup, int *ilp, int *iup, int *mp, 
	    double *W)
{
  int    info, itmp, j;
  bool   onlyW  = (jobz[0]  == 'N' || jobz[0]  == 'n');
  bool   wantZ  = (jobz[0]  == 'V' || jobz[0]  == 'v');
  bool   cntval = (jobz[0]  == 'C' || jobz[0]  == 'c');
  bool   upper  = (uplo[0]  == 'U' || uplo[0]  == 'u');
  bool   lower  = (uplo[0]  == 'L' || uplo[0]  == 'l');
  bool   alleig = (range[0] == 'A' || range[0] == 'a');
  bool   valeig = (range[0] == 'V' || range[0] == 'v');
  bool   indeig = (range[0] == 'I' || range[0] == 'i');
  int    n      = *np;
  int    lda    = *ldap;
  char   *trans;
  double complex *Z;
  double complex one = 1.0 + 0.0*I;

  /* Check input */
  if (n <= 1) return(1);
  if (*itype < 1 || *itype > 3) return(1);
  if (!lower  && !upper) return(1);
  if (!alleig && !valeig && !indeig) return(1);
  if (!onlyW  && !wantZ  && !cntval) return(1);

  if (indeig) itmp = *iup-*ilp+1;
  else        itmp = n;
  Z = (double complex*) malloc((size_t) itmp*n*sizeof(double complex));
  assert(Z != NULL);

  /* Form the Cholesky factor of B */
  zpotrf_(uplo, np, B, ldbp, &info);
  assert(info == 0);

  /* Convert problem to standard eigenvalue problem */
  zhegst_(itype, uplo, np, A, ldap, B, ldbp, &info);
  assert(info == 0);

  /* Solve standard eigenvalue problem using MRRR */
  info = zheeig(jobz, range, uplo, np, A, ldap, vlp, vup, 
		ilp, iup, mp, W, Z, np);
  assert(info == 0);

  /* Backtransform eigenvectors */
  if (wantZ) {

    if (*itype ==  1 || *itype == 2) {
      /* A*x = lambda*B*x or A*B*x = lambda*x require 
       * x = inv(L)'*y or x = inv(U)*y */

      if (upper) trans = "N";
      else       trans = "C";

      ztrsm_("Left", uplo, trans, "Non-unit", np, mp, &one, 
	     B, ldbp, Z, np);
    } else if (*itype == 3) {
      /* B*A*x = lambda*x requires x = L*y or U'*y */
      
      if (upper) trans = "C";
      else       trans = "N";

      ztrmm_("Left", uplo, trans, "Non-unit", np, mp, &one, 
	     B, ldbp, Z, np);
    } else {
      return(1);
    }

  } /* Backtransform eigenvectors */ 

  /* Copy eigenvectors into A */
  for (j=0; j<*mp; j++)
    memcpy(&A[j*lda], &Z[j*n], n*sizeof(double complex));

  free(Z);

  return(0);
}




/* Routine for the dense symmetric eigenproblem in packed storage */
int zhpeig(char *jobz, char *range, char *uplo, int *np, 
	   double complex *AP, double *vlp, double *vup, 
	   int *ilp, int *iup, int *mp, double *W, 
	   double complex *Z, int *ldzp)
{
  int      n      = *np;
  long int nn     = n;
  int      ldz    = *ldzp;
  int      tryRAC = 1;
  long int tmp, mm;
  long int i, j;

  double *D, *E, *Ztmp;
  double complex *TAU;
  double complex *work;
  int    *Zsupp;
  int    info, ione=1, m;

  double scale = 1.0;
  double invscale;

  bool   onlyW  = (jobz[0]  == 'N' || jobz[0]  == 'n');
  bool   wantZ  = (jobz[0]  == 'V' || jobz[0]  == 'v');
  bool   cntval = (jobz[0]  == 'C' || jobz[0]  == 'c');
  bool   alleig = (range[0] == 'A' || range[0] == 'a');
  bool   valeig = (range[0] == 'V' || range[0] == 'v');
  bool   indeig = (range[0] == 'I' || range[0] == 'i');

  if( !(onlyW  || wantZ  || cntval) ) return(1);
  if( !(alleig || valeig || indeig) ) return(1);
  if(n <= 0) return(1);
  if (valeig) {
    if(*vup<=*vlp) return(1);
  } else if (indeig) {
    if (*ilp<1 || *ilp>n || *iup<*ilp || *iup>n) return(1);
  }
  
  D = (double *) malloc( n*sizeof(double) );
  assert(D != NULL);
  
  E = (double *) malloc( n*sizeof(double) );
  assert(E != NULL);
  
  TAU = (double complex*) malloc( n*sizeof(double complex) );
  assert(TAU != NULL);
  
  work = (double complex*) malloc( n*sizeof(double complex) );
  assert(work != NULL);
  
  Zsupp = (int *) malloc( 2*n*sizeof(int) );
  assert(Zsupp != NULL);
  
  /* Scale matrix if necessary */
  scale = zscale_matrix_packed(range, uplo, np, AP, vlp, vup, work);

  /* Reduction to tridiagonal */
  zhptrd_(uplo, np, AP, D, E, TAU, &info);
  assert(info == 0);

  /* Use MRRR to compute eigenvalues and -vectors using part of Z 
   * to temporarily store the real eigenvectors of the tridiagonal */
  if (alleig)
    mm    = n;
  else if (indeig)
    mm    = (*iup)-(*ilp)+1;
  else {
    info = mrrr("Count", range, np, D, E, vlp, vup, ilp, iup, 
		&tryRAC, &m, W, NULL, ldzp, Zsupp);
    mm = m;
  }
  tmp  = (nn*mm)/2 + ( ((nn*mm) % 2) > 0 ); /* ceil(n*m/2) */
  Ztmp = (double *) &Z[ldz*mm - tmp];

  /* Actual call to MRRR */
  info = mrrr(jobz, range, np, D, E, vlp, vup, ilp, iup, 
  	      &tryRAC, mp, W, Ztmp, np, Zsupp);
  assert(info == 0);

  /* Copy intermediate real eigenvectors to complex Z */ 
  for (j=0; j<(*mp); j++)
    for (i=0; i<(*np); i++)
      Z[j*ldz + i] = Ztmp[j*n + i] + 0.0*I;
  
  /* Backtransformation Z = Q*Z */
  zupmtr_("L", uplo, "N", np, mp, AP, TAU, Z, ldzp, work, &info);
  assert(info == 0);

  /* Scaling of eigenvalues if necessary */
  if (scale != 1.0) { /* FP cmp okay */
    invscale = 1.0/scale;
    *vlp *= invscale;
    *vup *= invscale;
    dscal_(mp, &invscale, W, &ione);
  }
  
  free(D);
  free(E);
  free(TAU);
  free(work);
  free(Zsupp);
  
  return(0);
}




/* Routine for the dense generalized Hermitian-definite eigenproblem 
 * in packed storage */
int zhpgeig(int *itype, char *jobz, char *range, char *uplo, int *np, 
	    double complex *AP, double complex *BP, double *vlp, 
	    double *vup, int *ilp, int *iup, int *mp, double *W, 
	    double complex *Z, int *ldzp)
{
  int    n      = *np;
  int    ldz    = *ldzp;
  int    info, j, ione=1;
  bool   onlyW  = (jobz[0]  == 'N' || jobz[0]  == 'n');
  bool   wantZ  = (jobz[0]  == 'V' || jobz[0]  == 'v');
  bool   cntval = (jobz[0]  == 'C' || jobz[0]  == 'c');
  bool   upper  = (uplo[0]  == 'U' || uplo[0]  == 'u');
  bool   lower  = (uplo[0]  == 'L' || uplo[0]  == 'l');
  bool   alleig = (range[0] == 'A' || range[0] == 'a');
  bool   valeig = (range[0] == 'V' || range[0] == 'v');
  bool   indeig = (range[0] == 'I' || range[0] == 'i');
  char   *trans;

  /* Check input */
  if (n <= 1) return(1);
  if (*itype < 1 || *itype > 3) return(1);
  if (!lower  && !upper) return(1);
  if (!alleig && !valeig && !indeig) return(1);
  if (!onlyW  && !wantZ  && !cntval) return(1);

  /* Form the Cholesky factor of B */
  zpptrf_(uplo, np, BP, &info);
  assert(info == 0);

  /* Convert problem to standard eigenvalue problem */
  zhpgst_(itype, uplo, np, AP, BP, &info);
  assert(info == 0);

  /* Solve standard eigenvalue problem using MRRR */
  info = zhpeig(jobz, range, uplo, np, AP, vlp, vup, 
  		ilp, iup, mp, W, Z, ldzp);
  assert(info == 0);

  /* Backtransform eigenvectors */
  if (wantZ) {

    if (*itype ==  1 || *itype == 2) {
      /* A*x = lambda*B*x or A*B*x = lambda*x require 
       * x = inv(L)'*y or x = inv(U)*y */

      if (upper) trans = "N";
      else       trans = "C";

      for (j=0; j<*mp; j++)
        ztpsv_(uplo, trans, "Non-unit", np, BP, &Z[j*ldz], &ione); 

    } else if (*itype == 3) {
      /* B*A*x = lambda*x requires x = L*y or U'*y */
      
      if (upper) trans = "C";
      else       trans = "N";

      for (j=0; j<*mp; j++)
        ztpmv_(uplo, trans, "Non-unit", np, BP, &Z[j*ldz], &ione);
      
    } else {
      return(1);
    }

  } /* Backtransform eigenvectors */ 

  return(0);
}



static 
double zscale_matrix(char *range, char *uplo, int *np, 
		     double complex *A, int *ldap, double *vlp, 
		     double *vup, double complex *work)
{
  double sigma = 1.0;
  double smlnum, bignum, rmin, rmax;
  double norm;
  bool   scaled = false;
  bool   valeig = (range[0] == 'V' || range[0] == 'v');
  bool   lower  = (uplo[0] == 'L' || uplo[0] == 'l');
  int    n      = *np;
  int    lda    = *ldap;
  int    i, ione=1, itmp;

  smlnum = DBL_MIN / DBL_EPSILON;
  bignum = 1.0 / smlnum;
  rmin   = sqrt(smlnum);
  rmax   = fmin(sqrt(bignum), 1.0 / sqrt(sqrt(DBL_MIN)));

  norm = zlanhe_("M", uplo, np, A, ldap, work);
  if (norm > 0.0 && norm < rmin) {
    scaled = true;
    sigma  = rmin / norm;
  } else if (norm > rmax) {
    scaled = true;
    sigma  = rmax / norm;
  }
  if (scaled) {
    if (lower) {
      for (i=0; i<n; i++) {
	itmp = n - i;
	zdscal_(&itmp, &sigma, &A[i + i*lda], &ione);
      }
    } else {
      for (i=0; i<n; i++) {
	itmp = i + 1;
	zdscal_(&itmp, &sigma, &A[i*lda], &ione);
      }
    }
    if (valeig) {
      *vlp *= sigma;
      *vup *= sigma;
    }
  }
  
  return(sigma);
}



/* Scale matrix in packed storage to allowable range */ 
static 
double zscale_matrix_packed(char *range, char *uplo, int *np, 
			    double complex *AP, double *vlp, 
			    double *vup, double complex *work)
{
  double sigma = 1.0;
  double smlnum, bignum, rmin, rmax;
  double norm;
  bool   scaled = false;
  bool   valeig = (range[0] == 'V' || range[0] == 'v');
  int    n      = *np;
  int    ione=1, itmp;
  double *tmpp;

  smlnum = DBL_MIN / DBL_EPSILON;
  bignum = 1.0 / smlnum;
  rmin   = sqrt(smlnum);
  rmax   = fmin(sqrt(bignum), 1.0 / sqrt(sqrt(DBL_MIN)));

  tmpp = (double*) work;
  norm = zlanhp_("M", uplo, np, AP, tmpp);
  if (norm > 0.0 && norm < rmin) {
    scaled = true;
    sigma  = rmin / norm;
  } else if (norm > rmax) {
    scaled = true;
    sigma  = rmax / norm;
  }
  if (scaled) {
    itmp = (n*(n+1))/2;
    zdscal_(&itmp, &sigma, AP, &ione);
    if (valeig) {
      *vlp *= sigma;
      *vup *= sigma;
    }
  }
  
  return(sigma);
}
#endif




/* Fortran prototypes */
void dsyeig_(char *jobz, char *range, char *uplo, int *n, double *A, 
	     int *lda, double *vl, double *vu, int *il, int *iu, 
	     int *m, double *W, double *Z, int *ldz, int *info)
{
  *info = dsyeig(jobz, range, uplo, n, A, lda, vl, vu, il, iu, 
		 m, W, Z, ldz);
}


void dsygeig_(int *itype, char *jobz, char *range, char *uplo, 
	      int *np, double *A, int *ldap, double *B, int *ldbp, 
	      double *vlp, double *vup, int *ilp, int *iup, int *mp, 
	      double *W, int *info)
{
  *info = dsygeig(itype, jobz, range, uplo, np, A, ldap, B, ldbp, 
		  vlp, vup, ilp, iup, mp, W);
}


void dspeig_(char *jobz, char *range, char *uplo, int *np, double *AP, 
	     double *vlp, double *vup, int *ilp, int *iup, 
	     int *mp, double *W, double *Z, int *ldzp, int *info)
{
  *info = dspeig(jobz, range, uplo, np, AP, vlp, vup, ilp, iup, 
		 mp, W, Z, ldzp);
}


void dspgeig_(int *itype, char *jobz, char *range, char *uplo, int *np, 
	      double *AP, double *BP, double *vlp, double *vup, 
	      int *ilp, int *iup, int *mp, double *W, double *Z, 
	      int *ldzp, int *info)
{
  *info =  dspgeig(itype, jobz, range, uplo, np, AP, BP, vlp, vup, 
		   ilp, iup, mp, W, Z, ldzp);
}


#ifdef COMPLEX_SUPPORTED
void zheeig_(char *jobz, char *range, char *uplo, int *n, 
	     double complex *A, int *lda, double *vl, double *vu, 
	     int *il, int *iu, int *m, double *W, double complex *Z, 
	     int *ldz, int *info)
{
  *info = zheeig(jobz, range, uplo, n, A, lda, vl, vu, il, iu, 
		 m, W, Z, ldz);
}


void zhegeig_(int *itype, char *jobz, char *range, char *uplo, int *np, 
	      double complex *A, int *ldap, double complex *B, int *ldbp, 
	      double *vlp, double *vup, int *ilp, int *iup, int *mp, 
	      double *W, int *info)
{
  *info = zhegeig(itype, jobz, range, uplo, np, A, ldap, B, ldbp, 
		  vlp, vup, ilp, iup, mp, W);
}


void zhpeig_(char *jobz, char *range, char *uplo, int *np, 
	     double complex *AP, double *vlp, double *vup, 
	     int *ilp, int *iup, int *mp, double *W, 
	     double complex *Z, int *ldzp, int *info)
{
  *info =  zhpeig(jobz, range, uplo, np, AP, vlp, vup, 
		  ilp, iup, mp, W, Z, ldzp);
}


void zhpgeig_(int *itype, char *jobz, char *range, char *uplo, int *np, 
	      double complex *AP, double complex *BP, double *vlp, 
	      double *vup, int *ilp, int *iup, int *mp, double *W, 
	      double complex *Z, int *ldzp, int *info)
{
  *info =  zhpgeig(itype, jobz, range, uplo, np, AP, BP, vlp, 
		   vup, ilp, iup, mp, W, Z, ldzp);
}
#endif

