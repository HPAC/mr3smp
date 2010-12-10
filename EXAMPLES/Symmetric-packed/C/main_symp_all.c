/* Example code for how to call the routine "dspeig" for the dense 
 * symmetric eigenproblem with matrix A in packed storage.
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
#include <assert.h>
#include "global.h"
#include "mrrr.h"


int main(int argc, char **argv)
{
  int n   = 100;
  int ldz = n;
  int nz  = n;
  int size = (n*(n+1))/2;

  double *AP, *W, *Z;
  double vl, vu;
  int    il, iu;
  int    m, i;

  AP = (double *) malloc((size_t) size*sizeof(double));
  assert(AP != NULL);

  W = (double *) malloc((size_t)   n *sizeof(double));
  assert(W != NULL);

  Z = (double *) malloc((size_t) n*nz*sizeof(double));
  assert(Z != NULL);

  /* Initialize matrix in packed storage; since we use random value 
   * entries no distinction between upper or lower part stored */  
  srand( (unsigned) time( (time_t *) NULL) );
  for (i=0; i<size; i++)
    AP[i] = (double) rand() / RAND_MAX;


  /* Calling "dspeig" to compute all or a subset of eigenvalues and 
   * optinally eigenvectors of a dense symmetric matrix in packed storage 
   * using LAPACK routines for reduction and backtransformation and 
   * the multi-threaded MRRR for the tridiagonal stage.
   * The number of threads for the LAPACK routines are set by 
   * OMP_NUM_THREADS or GOTO_NUM_THREADS or MKL_NUM_THREADS or ... 
   * depending on the BLAS used. For the tridiagonal stage with 
   * PMR_NUM_THREADS. 
   * Since A is initialized randomly we can either assume it represents 
   * the upper or lower triangular part of A in packed storage */
  dspeig("Vectors", "All", "Lower", &n, AP, &vl, &vu, &il, &iu,
	 &m, W, Z, &ldz);


  printf("Successfully computed eigenpairs!\n");

  free(AP);
  free(W);
  free(Z);

  return(0);
}
