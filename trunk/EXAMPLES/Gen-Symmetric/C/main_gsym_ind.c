/* Example code for how to call the routine "dsygeig" for the dense 
 * symmetric-definite generalized eigenproblem.
 * The routine must be linked to LAPACK and an optimized BLAS.
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

static void init_symmetric_matrix(double*, int, int);



int main(int argc, char **argv)
{
  int    i, info;
  int    n   = 100;
  int    lda = n;
  int    ldb = n;
  double *A, *B, *W;
  double vl, vu;
  int    il, iu;
  int    m;
  int    itype;

  A = (double *) malloc((size_t) n*n *sizeof(double));
  assert(A != NULL);

  B = (double *) malloc((size_t) n*n *sizeof(double));
  assert(B != NULL);

  W = (double *) malloc((size_t)   n *sizeof(double));
  assert(W != NULL);

  /* Initialize symmetric matrix fully, that is upper and lower 
   * triangular part are inbitialized */
  init_symmetric_matrix(A, n, lda);
  init_symmetric_matrix(B, n, ldb);

  /* Making B positive definite */
  for (i=0; i<n; i++)
    B[i+i*ldb] += n; 


  /* Calling "dsygeig" to compute all or a subset of eigenvalues and 
   * optinally eigenvectors of a dense generalized symmetric-definite 
   * eigenproblem;
   * The number of threads for the LAPACK routines are set by 
   * OMP_NUM_THREADS or GOTO_NUM_THREADS or MKL_NUM_THREADS or ... 
   * depending on the BLAS used. For the tridiagonal stage with 
   * PMR_NUM_THREADS. 
   * Since A and B are initialized fully, both the upper or lower 
   * triangular parts can be used. */
  itype = 1; /* A*x = lambda*B*x */
  il    = 1; /* To compute eigenpairs 'il' to 'iu' */ 
  iu    = n/2 + 1;
  info = dsygeig(&itype, "Vectors", "Index", "Low", &n, A, &lda, 
		 B, &ldb, &vl, &vu, &il, &iu, &m, W);
  assert(info == 0);


  printf("Successfully computed eigenpairs!\n");

  free(A);
  free(B);
  free(W);

  return(0);
}




/* Auxiliary routine to initialize matrix */
static void init_symmetric_matrix(double *A, int n, int lda)
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
