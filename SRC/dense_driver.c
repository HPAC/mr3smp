/* Routines for the dense symmetric and hermitian eigenproblem 
 * using the multithreaded MRRR for the tridiagonal stage. */

#include <stdlib.h>
#include <assert.h>
#include <complex.h>
#include "global.h"
#include "mrrr.h"  


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
  int    info;

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

  TAU = (double *) malloc( n*n*sizeof(double) );
  assert(TAU != NULL);

  work = (double *) malloc( szwrk*sizeof(double) );
  assert(work != NULL);

  Zsupp = (int *) malloc( 2*n*sizeof(int) );
  assert(Zsupp != NULL);


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


  free(D);
  free(E);
  free(TAU);
  free(work);
  free(Zsupp);

  return(0);
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

  TAU = (double complex *) malloc((size_t) n*n*sizeof(double complex));
  assert(TAU != NULL);

  work = (double complex *) malloc( szwrk*sizeof(double complex) );
  assert(work != NULL);

  Zsupp = (int *) malloc(2*n*sizeof(int));
  assert(Zsupp != NULL);


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
  for (i=0; i<(*np); i++)
    for (j=0; j<(*mp); j++)
      Z[ i * (*ldzp) + j ] = Ztmp[ i * (*np) + j ] + 0.0 * I;
  
  /* Backtransformation Z = U*Z */
  zunmtr_("L", uplo, "N", np, mp, A, ldap, TAU, Z, ldzp, 
	  work, &szwrk, &info);
  assert(info == 0);


  free(D);
  free(E);
  free(TAU);
  free(work);
  free(Zsupp);

  return(0);
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



#ifdef COMPLEX_SUPPORTED
void zheeig_(char *jobz, char *range, char *uplo, int *n, 
	     double complex *A, int *lda, double *vl, double *vu, 
	     int *il, int *iu, int *m, double *W, double complex *Z, 
	     int *ldz, int *info)
{
  *info = zheeig(jobz, range, uplo, n, A, lda, vl, vu, il, iu, 
		 m, W, Z, ldz);
}
#endif

