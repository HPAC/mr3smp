/* Example file for using MRRR to compute a subset of eigenpairs 
 * specified by value.
 *
 * Example run: 
 * $ export OMP_NUM_THREADS=4; ./main_val.x
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "mrrr.h"

static int  read_tri_mat(char*, double**, double**);
static void print_vector(char*, double*, char*, int);
static void print_matrix(char*, double*, int, int, int);


int main(int argc, char **argv)
{
  /* Input parameter to 'mrrr' */
  int     n;              /* Matrix size */
  int     il, iu;         
  int     tryRAC = 1;     /* Try high rel. accuracy */
  double  *D, *E;         /* Diagonal and off-diagonal elements */
  double  vl, vu;         

  double  *W;             /* eigenvalues */
  int     m;              /* number of eigenpairs */
  double  *Z;             /* eigenvectors; stored by colums */
  int     ldz;
  int     *Zsupp;         /* eigenvector support */

  /* Others */
  int     info;

  /* Read in data from file, space for D and E will
   * be allocated and needs to be freed at the end */
  n = read_tri_mat("./Wilkinson21.data", &D, &E);
  ldz = n;

  /* Print input */
  printf("\n%% Input matrix:\n\n");
  printf("n = %d;\n", n);
  print_vector("D=[", D, "];", n  );
  print_vector("E=[", E, "];", n-1);

  /* Compute eigenpairs contained in (vl,vu] */
  vl = -1.22;
  vu =  5.93;

  /* Allocate memory */
  W     = (double *) malloc( n    * sizeof(double) );
  Zsupp = (int *)    malloc( 2*n  * sizeof(int)    );

  info = mrrr("Count", "Value", &n, D, E, &vl, &vu, &il, &iu,
	      &tryRAC, &m, W, NULL, &ldz, Zsupp);
  assert(info == 0);
  
  Z     = (double *) malloc((size_t) n*m * sizeof(double) );

  /* Use MRRR to compute eigenvalues and -vectors */
  info = mrrr("Vectors", "Value", &n, D, E, &vl, &vu, &il,
	      &iu, &tryRAC, &m, W, Z, &ldz, Zsupp);
  assert(info == 0);

  /* Print results */
  printf("\n\n%% Results:\n");
  print_vector("W=[", W, "];", m);
  print_matrix("Z", Z, n, m, 0);

  /* Free allocated memory */
  free(D);
  free(E);
  free(W);
  free(Z);
  free(Zsupp);

  return(0);
}




/* 
 * Reads the triadiagonal matrix from a file.
 */
static int read_tri_mat(char *filename, double **Dp, double **Ep)
{
  int    i, n;
  FILE   *filedes;

  filedes = fopen(filename, "r");
  assert(filedes != NULL);

  fscanf(filedes, "%d", &n);

  *Dp = (double *) malloc( n * sizeof(double) );
  assert(*Dp != NULL);

  *Ep = (double *) malloc( n * sizeof(double) );
  assert(*Ep != NULL);

  for (i=0; i<n; i++) {
    fscanf(filedes, "%le %le", *Dp+i, *Ep+i);
  }
  (*Ep)[n-1] = 0.0;

  fclose(filedes);

  return(n);
}



static void print_vector(char *pre, double *v, char *post, int n)
{
  int i;

  printf("\n%s\n", pre);
  for (i=0; i<n; i++) {
    printf("%.17e", v[i]);
    printf("\n");
  }
  printf("%s\n", post);
}



static void print_matrix(char *name, double *A, int lda, 
			 int numcols, int offset)
{
  int i, j;
  char str1[20], str2[20];

  for (j=1; j<=numcols; j++) {
    str1[0] = '\0';
    strcat(str1,name);
    strcat(str1,"(:,");
    i = sprintf(str2, "%d",j+offset);
    strcat(str1, str2);
    strcat(str1, ")=[");
    print_vector(str1, &A[(j-1)*lda],"];", lda);
  }
}
