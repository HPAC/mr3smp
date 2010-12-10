/* Example code for how to call the routine "zhpeig" for the dense 
 * Hermitian eigenproblem with matrix A in packed storage.
 * The routine uses LAPACK routines for the reduction and 
 * back-transformation and must therefore be linked to LAPACK.
 * The number of threads for the LAPACK routines are set by 
 * OMP_NUM_THREADS or GOTO_NUM_THREADS or MKL_NUM_THREADS or ... 
 * depending on the BLAS used. For the tridiagonal stage with 
 * PMR_NUM_THREADS. 
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <complex.h> 
#include <assert.h>
#include "global.h"
#include "mrrr.h"



int main(int argc, char **argv)
{
  int n       = 100;
  int ldz     = n;
  int nz      = n;
  size_t size = (n*(n+1))/2;

  double complex *AP, *Z;
  double         *W;
  double         vl, vu;
  int            il, iu;
  int            m, i, j, k;
  double         tmp1, tmp2;

  AP = (double complex*) malloc(size*sizeof(double complex));
  assert(AP != NULL);

  W = (double *) malloc((size_t)   n *sizeof(double));
  assert(W != NULL);

  Z = (double complex*) malloc((size_t) n*nz*sizeof(double complex));
  assert(Z != NULL);

  /* Initialize matrix in packed storage; assume upper part stored */  
  srand((unsigned) time( (time_t *) NULL));
  k = 0;
  for (j=0; j<n; j++) {
    for (i=0; i<=j; i++) {
      tmp1 = (double) rand() / RAND_MAX;
      tmp2 = (double) rand() / RAND_MAX;
      if (i==j) AP[k++] = tmp1 + 0.0*I;
      else AP[k++] = tmp1 + tmp2*I;
    }
  }



  /* Calling "zhpeig" to compute all or a subset of eigenvalues and 
   * optinally eigenvectors of a dense Hermitian matrix in packed storage 
   * using LAPACK routines for reduction and backtransformation and 
   * the multi-threaded MRRR for the tridiagonal stage.
   * The number of threads for the LAPACK routines are set by 
   * OMP_NUM_THREADS or GOTO_NUM_THREADS or MKL_NUM_THREADS or ... 
   * depending on the BLAS used. For the tridiagonal stage with 
   * PMR_NUM_THREADS. 
   * Since A is assumed to store the upper triangular part of the
   * matrix we must call the routine with UPLO="U" */
  il = 1;
  iu = n/2 + 1;
  zhpeig("Vectors", "Index", "Upper", &n, AP, &vl, &vu, &il, &iu,
	 &m, W, Z, &ldz);
  
  
  printf("Successfully computed eigenpairs!\n");

  free(AP);
  free(W);
  free(Z);

  return(0);
}
