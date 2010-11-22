/* Example code for how to call the routine "dsyeig" for the dense 
 * symmetric eigenproblem.
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

void init_symmetric_matrix(double*, int, int);



int main(int argc, char **argv)
{
  int n   = 300;
  int lda = n;
  int ldz = n;
  int nz  = n;

  double *A, *W, *Z;
  double vl, vu;
  int    il, iu;
  int    m;

  A = (double *) malloc((size_t) n*n *sizeof(double));
  assert(A != NULL);

  W = (double *) malloc((size_t)   n *sizeof(double));
  assert(W != NULL);

  Z = (double *) malloc((size_t) n*nz*sizeof(double));
  assert(Z != NULL);


  /* Initialize symmetric matrix fully, that is upper and lower 
   * triangular part are inbitialized */
  init_symmetric_matrix(A,n,lda);

  /* Calling "dsyeig" to compute all or a subset of eigenvalues and 
   * optinally eigenvectors of a dense symmetric matrix using LAPACK 
   * routines for reduction and backtransformation and the multi-threaded 
   * MRRR for the tridiagonal stage.
   * The number of threads for the LAPACK routines are set by 
   * OMP_NUM_THREADS or GOTO_NUM_THREADS or MKL_NUM_THREADS or ... 
   * depending on the BLAS used. For the tridiagonal stage with 
   * PMR_NUM_THREADS. 
   * Since A is initialized fully, both the upper or lower triangular part 
   * of A can be used in the call to "dsyeig". */
  dsyeig("Vectors", "All", "Low", &n, A, &lda, &vl, &vu, &il, &iu, 
	 &m, W, Z, &ldz);


  printf("Successfully computed eigenpairs!\n");

  free(A);
  free(W);
  free(Z);

  return(0);
}




/* Auxiliary routine to initialize matrix */
void init_symmetric_matrix(double *A, int n, int lda)
{
  int i, j;

  srand( (unsigned) time( (time_t *) NULL) );

  for (j=0; j<n; j++)
    for (i=j; i<n; i++)
      A[ i + j*lda ] = (double) rand() / RAND_MAX;

  for (j=0; j<n; j++)
    for (i=0; i<j; i++)
      A[ i + j*lda ] = A[ j + i*lda ];
}
