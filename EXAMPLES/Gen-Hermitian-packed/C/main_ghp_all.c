/* Example code for how to call the routine "zhpgeig" for the dense 
 * hermitian-definite generalized eigenproblem in packed storage format.
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
#include <complex.h>
#include "global.h"
#include "mrrr.h"

static void init_herm_matrix_packed(double complex*, int);



int main(int argc, char **argv)
{
  int            i, j, k, itype, info;
  int            n    = 100;
  int            size = (n*(n+1))/2;
  int            ldz  = n;

  double complex *AP, *BP, *Z;
  double         *W;
  double         vl, vu;
  int            il, iu;
  int            m;

  AP = (double complex *) malloc((size_t) size*sizeof(double complex));
  assert(AP != NULL);

  BP = (double complex *) malloc((size_t) size*sizeof(double complex));
  assert(BP != NULL);

  W = (double *) malloc((size_t)   n *sizeof(double));
  assert(W != NULL);

  Z = (double complex *) malloc((size_t) n*n*sizeof(double complex));
  assert(Z != NULL);

  /* Initialize hermitian matrix fully, that is upper and lower 
   * triangular part are inbitialized */
  init_herm_matrix_packed(AP, n);
  init_herm_matrix_packed(BP, n);

  /* Making B positive definite, assuming storage of upper part */
  k = 0;
  for (j=0; j<n; j++)
    for (i=0; i<=j; i++) {
      if (i==j) BP[k] += n;
      k++;
    }



  /* Calling "zhpgeig" to compute all or a subset of eigenvalues and 
   * optinally eigenvectors of a dense generalized symmetric-definite 
   * eigenproblem;
   * The number of threads for the LAPACK routines are set by 
   * OMP_NUM_THREADS or GOTO_NUM_THREADS or MKL_NUM_THREADS or ... 
   * depending on the BLAS used. For the tridiagonal stage with 
   * PMR_NUM_THREADS. 
   * Since A and B are initialized fully, both the upper or lower 
   * triangular parts can be used. */
  itype = 1; /* A*x = lambda*B*x */
  info = zhpgeig(&itype, "Vectors", "All", "Upper", &n, AP, 
		 BP, &vl, &vu, &il, &iu, &m, W, Z, &ldz);
  assert(info == 0);



  printf("Successfully computed eigenpairs!\n");

  free(AP);
  free(BP);
  free(Z);
  free(W);

  return(0);
}




/* Auxiliary routine to initialize matrix, such that the diagonal
 * is real for cases of upper and lower part stored */
static void init_herm_matrix_packed(double complex *AP, int n)
{
  int i, j, k;
  double tmp1, tmp2;

  srand( (unsigned) time( (time_t *) NULL) );

  k = 0;
  for (j=0; j<n; j++) {
    for (i=j; i<n; i++) {
      tmp1 = (double) rand() / RAND_MAX;
      tmp2 = (double) rand() / RAND_MAX;
      if (i==j) AP[k++] = tmp1 + 0.0*I;
      else      AP[k++] = tmp1 + tmp2*I;
    }
  }

  k = 0;
  for (j=0; j<n; j++)
    for (i=0; i<=j; i++) {
      if (i==j) AP[k] = creal(AP[k]) + 0.0*I;
      k++;
    }
}
