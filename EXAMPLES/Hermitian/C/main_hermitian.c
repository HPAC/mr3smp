/* Example code for how to call the routine "dsyeig" for the dense 
 * hermitian eigenproblem.
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
#include <complex.h>
#include "global.h"
#include "mrrr.h"

static void init_hermitian_matrix(double complex*, int, int);



int main(int argc, char **argv)
{
  int            n   = 100;
  int            lda = n;
  int            ldz = n;
  int            nz  = n;

  double complex *A, *Z;
  double         *W;
  double         vl, vu;
  int            il, iu;
  int            m;

  A = (double complex *) malloc((size_t) n*n *sizeof(double complex));
  assert(A != NULL);

  W = (double *) malloc((size_t)   n *sizeof(double));
  assert(W != NULL);

  Z = (double complex *) malloc((size_t) n*nz*sizeof(double complex));
  assert(Z != NULL);

  /* Initialize hermitian matrix fully, that is upper and lower 
   * triangular part are inbitialized */
  init_hermitian_matrix(A,n,lda);


  /* Calling "dsyeig" to compute all or a subset of eigenvalues and 
   * optinally eigenvectors of a dense hermitian matrix using LAPACK 
   * routines for reduction and backtransformation and the multi-threaded 
   * MRRR for the tridiagonal stage.
   * The number of threads for the LAPACK routines are set by 
   * OMP_NUM_THREADS or GOTO_NUM_THREADS or MKL_NUM_THREADS or ... 
   * depending on the BLAS used. For the tridiagonal stage with 
   * PMR_NUM_THREADS. 
   * Since A is initialized fully, both the upper or lower triangular part 
   * of A can be used in the call to "dsyeig". */
  zheeig("Vectors", "All", "Low", &n, A, &lda, &vl, &vu, &il, &iu, 
	 &m, W, Z, &ldz);


  printf("Successfully computed eigenpairs!\n");

  free(A);
  free(W);
  free(Z);

  return(0);
}




/* Auxiliary routine to initialize matrix */
static void init_hermitian_matrix(double complex *A, int n, int lda)
{
  int i, j;
  double tmp1, tmp2;

  srand( (unsigned) time( (time_t *) NULL) );

  for (j=0; j<n; j++) {
    for (i=j; i<n; i++) {
      tmp1 = (double) rand() / RAND_MAX;
      tmp2 = (double) rand() / RAND_MAX;
      A[ i + j*lda ] = tmp1 + tmp2*I;
    }
  }

  for (j=0; j<n; j++) {
    for (i=0; i<j; i++)
      A[ i + j*lda ] = conj(A[ j + i*lda ]);
    A[ j + j*lda ] = creal(A[ j + j*lda ]) + 0.0*I;
  }
}
