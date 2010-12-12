/* Example code for how to call the routine "dspgeig" for the dense 
 * symmetric-definite generalized eigenproblem in packed storage format.
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

static void init_sym_matrix_packed(double*, int);



int main(int argc, char **argv)
{
  int    i, j ,k, info;
  int    n    = 100;
  int    size = (n*(n+1))/2;
  double *AP, *BP, *W, *Z;
  double vl, vu;
  int    il, iu;
  int    m;
  int    itype;
  int    ldz = n;

  AP = (double *) malloc((size_t) size*sizeof(double));
  assert(AP != NULL);

  BP = (double *) malloc((size_t) size*sizeof(double));
  assert(BP != NULL);

  W = (double *) malloc((size_t)    n *sizeof(double));
  assert(W != NULL);

  Z = (double *) malloc((size_t)  n*n *sizeof(double));
  assert(Z != NULL);

  /* Initialize symmetric matrix fully, that is upper and lower 
   * triangular part are inbitialized */
  init_sym_matrix_packed(AP, n);
  init_sym_matrix_packed(BP, n);

  /* Making B positive definite, assuming storage of upper part */
  k = 0;
  for (j=0; j<n; j++)
    for (i=0; i<=j; i++) {
      if (i==j) BP[k] += n;
      k++;
    }


  /* Calling "dspgeig" to compute all or a subset of eigenvalues and 
   * optinally eigenvectors of a dense generalized symmetric-definite 
   * eigenproblem;
   * The number of threads for the LAPACK routines are set by 
   * OMP_NUM_THREADS or GOTO_NUM_THREADS or MKL_NUM_THREADS or ... 
   * depending on the BLAS used. For the tridiagonal stage with 
   * PMR_NUM_THREADS. 
   * Since B is initialized assuming upper parts are stored, the 
   * routine is called with UPLO="U". */
  itype = 1; /* A*x = lambda*B*x */
  info = dspgeig(&itype, "Vectors", "All", "Upper", &n, AP,
		 BP, &vl, &vu, &il, &iu, &m, W, Z, &ldz);
  assert(info == 0);


  printf("Successfully computed eigenpairs!\n");

  free(AP);
  free(BP);
  free(Z);
  free(W);

  return(0);
}




/* Auxiliary routine to initialize matrix */
static void init_sym_matrix_packed(double *AP, int n)
{
  int i, j, k;

  srand( (unsigned) time( (time_t *) NULL) );

  k = 0;
  for (j=0; j<n; j++)
    for (i=j; i<n; i++)
      AP[k++] = (double) rand() / RAND_MAX;
}
