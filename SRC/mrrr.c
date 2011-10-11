/* Computing all or a subset of eigenvalues and optionally eigenvectors
 * of a symmetric tridiagonal matrix.
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
 * September 2010, Version 1.1
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
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "global.h"
#include "mrrr.h"
#include "structs.h"


int mrrr_val(char *jobz, char *range, int nthreads, in_t *Dstruct,
	     val_t *Wstruct, tol_t *tolstruct);
int mrrr_vec(int nthreads, in_t *Dstruct, val_t *Wstruct, 
	     vec_t *Zstruct, tol_t *tolstruct);

static inline void clean_up(double*, double*, double*, int*, int*, 
			    int*, int*, in_t*, val_t*, vec_t*, 
			    tol_t*);
static inline void scale_matrix(in_t*, val_t*, bool, double*);
static inline void invscale_eigvals(val_t*, double);
static int handle_small_cases(char*, char*, int*, double*, double*, 
			      double*, double*, int*, int*, int*, 
			      int*, double*, double*, int*, int*);
static void refine_to_highrac(double*, double*, in_t*, val_t*, 
			      tol_t*);
static inline void sort_eigenpairs(val_t *Wstruct, vec_t *Zstruct,
				   double *work);
static int cmpa(const void *, const void *);
static int cmpb(const void *, const void *);





int mrrr(char *jobz, char *range, int *np, double *restrict D,
	 double *restrict E, double *vlp, double *vup, int *ilp, 
	 int *iup, int *tryracp, int *mp, double *restrict W, 
	 double *restrict Z, int *ldz, int *restrict Zsupp)
{
  int    n      = *np;
  bool   onlyW  = (jobz[0]  == 'N' || jobz[0]  == 'n');
  bool   wantZ  = (jobz[0]  == 'V' || jobz[0]  == 'v');
  bool   cntval = (jobz[0]  == 'C' || jobz[0]  == 'c');
  bool   alleig = (range[0] == 'A' || range[0] == 'a');
  bool   valeig = (range[0] == 'V' || range[0] == 'v');
  bool   indeig = (range[0] == 'I' || range[0] == 'i');

  char   *ompvar;
  int    nthreads;

  double *Dcopy, *E2copy;
  double *gersch;
  double *Wgap, *Werr;
  int    *isplit, *iblock;
  int    *Windex, *Zindex;

  in_t   *Dstruct;  
  val_t  *Wstruct;
  vec_t  *Zstruct;
  tol_t  *tolstruct;

  double scale;              
  int    i, info;                    
  int    iblk, ishift;

  /* Check input parameters */
  if( !(onlyW  || wantZ  || cntval) ) return(1);
  if( !(alleig || valeig || indeig) ) return(1);
  if(n <= 0) return(1);
  if (valeig) {
    if(*vup<=*vlp) return(1);
  } else if (indeig) {
    if (*ilp<1 || *ilp>n || *iup<*ilp || *iup>n) return(1);
  }

  /* If only maximal number of local eigenvectors are queried
   * return if possible here */
  if (cntval) {
    if ( alleig ) {
      *mp = n;
      return(0);
    } else if (indeig) {
      *mp = *iup - *ilp + 1;
      return(0);
    }
  }

  /* Determine number of threads to be used */
  ompvar = getenv("PMR_NUM_THREADS");
  if (ompvar == NULL) {
    nthreads = DEFAULT_NUM_THREADS;
  } else {
    nthreads = atoi(ompvar);
  }

  /* Handle sequntial case and small cases by calling dstemr */
  if (n < DSTEMR_IF_SMALLER || nthreads == 1) {
    info = handle_small_cases(jobz, range, np, D, E, vlp, vup, ilp,
			      iup, tryracp, mp, W, Z, ldz, Zsupp);
    return(info);
  }

  /* Allocate memory shared by all threads */
  Werr   = (double *) malloc( n * sizeof(double) );
  assert(Werr != NULL);
  Wgap   = (double *) malloc( n * sizeof(double) );
  assert(Wgap != NULL);
  gersch = (double *) malloc( 2*n*sizeof(double) );
  assert(gersch != NULL);
  /* use calloc to initialize iblock, because for case range="I" and 
   * a block size of 1 the check if(range=I && iblock[index]==j) will 
   * use an uninitialized value otherwise */
  iblock = (int *)    calloc( n , sizeof( int  ) );
  assert(iblock != NULL);
  Windex = (int *)    malloc( n * sizeof( int  ) );
  assert(Windex != NULL);
  isplit = (int *)    malloc( n * sizeof( int  ) );
  assert(isplit != NULL);
  Zindex = (int *)    malloc( n * sizeof( int  ) );
  assert(Zindex != NULL);
  Dstruct   = (in_t *)   malloc( sizeof(in_t) );
  assert(Dstruct != NULL);
  Wstruct   = (val_t *)  malloc( sizeof(val_t) );
  assert(Wstruct != NULL);
  Zstruct   = (vec_t *)  malloc( sizeof(vec_t) );
  assert(Zstruct != NULL);
  tolstruct = (tol_t *)  malloc( sizeof(tol_t) );
  assert(tolstruct != NULL);

  /* Bundle variables into a structures */
  Dstruct->n               = n;
  Dstruct->D               = D;
  Dstruct->E               = E;
  Dstruct->isplit          = isplit;

  Wstruct->n               = n;
  Wstruct->vlp             = vlp;
  Wstruct->vup             = vup;
  Wstruct->ilp             = ilp;
  Wstruct->iup             = iup;
  Wstruct->mp              = mp;
  Wstruct->W               = W;
  Wstruct->Werr            = Werr;
  Wstruct->Wgap            = Wgap;
  Wstruct->Windex          = Windex;
  Wstruct->iblock          = iblock;
  Wstruct->gersch          = gersch;

  Zstruct->ldz             = *ldz;
  Zstruct->Z               = Z;
  Zstruct->Zsupp           = Zsupp;
  Zstruct->Zindex          = Zindex;

  /* Scale matrix if necessary */
  scale_matrix(Dstruct, Wstruct, valeig, &scale);

  /* Test if matrix warrants more expensive computations which
   * guarantees high relative accuracy */
  if (*tryracp) dlarrr_(&n, D, E, &info); /* 0 - rel acc */
  else info = -1;

  /* Set tolerance for splitting */ 
  if (info == 0) {
    tolstruct->split = DBL_EPSILON;
    /* Copy original diagonal, needed for refinement later */
    Dcopy = (double *) malloc( n * sizeof(double) );
    memcpy(Dcopy, D, n*sizeof(double));
    E2copy = (double *) malloc( n * sizeof(double) );
    assert(E2copy != NULL);
    for (i=0; i<n-1; i++) E2copy[i] = E[i]*E[i];
  } else {
    /* Use neg. threshold to force old splitting criterion */
    tolstruct->split = -DBL_EPSILON;
    *tryracp = 0;
  }

  /* Set tolerances */
  if (!wantZ) {
    tolstruct->rtol1 = 4.0 * DBL_EPSILON;
    tolstruct->rtol2 = 4.0 * DBL_EPSILON;
  } else {
    tolstruct->rtol1 = sqrt(DBL_EPSILON);
    tolstruct->rtol1 = fmin(tolstruct->rtol1, 1e-2 * MIN_RELGAP);
    tolstruct->rtol2 = sqrt(DBL_EPSILON)*5.0E-3;
    tolstruct->rtol2 = fmin(tolstruct->rtol2, MIN_RELGAP * 5.0E-6);
    tolstruct->rtol2 = fmax(tolstruct->rtol2, 4.0 * DBL_EPSILON);
  }
  /* LAPACK: tolstruct->bsrtol = sqrt(DBL_EPSILON); */
  tolstruct->bsrtol = fmin(tolstruct->rtol1, sqrt(DBL_EPSILON));
  tolstruct->RQtol  = 2.0 * DBL_EPSILON;
 
  /*  Compute the desired eigenvalues */
  info = mrrr_val(jobz, range, nthreads, Dstruct, Wstruct, tolstruct);
  assert(info == 0);

  /* If jobz="C" only the number of eigenpairs was computed */
  if (cntval) {    
    clean_up(Werr, Wgap, gersch, iblock, Windex, isplit, Zindex, 
	     Dstruct, Wstruct, Zstruct, tolstruct);
    return(0);
  }

  /*  Compute, if desired, associated eigenvectors */  
  if (wantZ) {

    info = mrrr_vec(nthreads, Dstruct, Wstruct, Zstruct, tolstruct);
    assert(info == 0);

  } else {
    /* Shift eigenvalues: 'mrrr_val' returns eigenvalues of SHIFTED 
     * root reprensentations; 'mrrr_vec' of UNshifted input T */
    for (i=0; i<*mp; i++) {
      iblk     = iblock[i];
      ishift   = isplit[iblk-1] - 1;
      W[i]    += E[ishift];           
    }
  } /* end of eigenvectors computed or not */

  /* Refine to high relative accuracy with respect to input T */
  if (*tryracp)
    refine_to_highrac(Dcopy, E2copy, Dstruct, Wstruct, tolstruct);
  
  /* If matrix was scaled, rescale eigenvalues */
  if (scale != 1.0) invscale_eigvals(Wstruct, scale); /* FP cmp OK */

  if (wantZ) {
    /* Sort eigenpairs; note: using 'gersch' as work space */
    sort_eigenpairs(Wstruct, Zstruct, gersch);
  } else {
    /* Sort eigenvalues only */
    qsort(W, *mp, sizeof(double), cmpa);
  }

  /* Free allocated memory */
  clean_up(Werr, Wgap, gersch, iblock, Windex, isplit, Zindex, 
	   Dstruct, Wstruct, Zstruct, tolstruct);
  if (*tryracp) {
    free(Dcopy);
    free(E2copy);
  }

  return(0);
}




/* Free's on allocated memory of pmrrr routine */
static inline   
void clean_up(double *Werr, double *Wgap, double *gersch, 
	      int *iblock, int *Windex, int *isplit, int *Zindex,  
	      in_t *Dstruct, val_t *Wstruct, vec_t *Zstruct,
	      tol_t *tolstruct)
{
  free(Werr);
  free(Wgap);
  free(gersch);
  free(iblock);
  free(Windex);
  free(isplit);
  free(Zindex);
  free(Dstruct);
  free(Wstruct);
  free(Zstruct);
  free(tolstruct);
}



/* Scale matrix to allowable range, returns 1.0 if not scaled */
static inline  
void scale_matrix(in_t *Dstruct, val_t *Wstruct, bool valeig, 
		  double *scale)
{
  int              n   = Dstruct->n;
  double *restrict D   = Dstruct->D;
  double *restrict E   = Dstruct->E;
  double          *vlp = Wstruct->vlp;
  double          *vup = Wstruct->vup;

  double           T_norm;              
  double           smlnum, bignum, rmin, rmax;
  int              IONE = 1, itmp;

  *scale = 1.0;

  /* set some machine dependent constants */
  smlnum = DBL_MIN / DBL_EPSILON;
  bignum = 1.0 / smlnum;
  rmin   = sqrt(smlnum);
  rmax   = fmin(sqrt(bignum), 1.0 / sqrt(sqrt(DBL_MIN)));

  T_norm = dlanst_("M", &n, D, E);  /* returns max(|T(i,j)|) */
  if (T_norm > 0 && T_norm < rmin) {
    *scale = rmin / T_norm;
  } else if (T_norm > rmax) {
    *scale = rmax / T_norm;
  }

  if (*scale != 1.0) {  /* FP comparison okay */
    /* Scale matrix and matrix norm */
    itmp = n-1;
    dscal_(&n,    scale, D, &IONE);
    dscal_(&itmp, scale, E, &IONE);

    if (valeig == true) {
      /* Scale eigenvalue bounds */
      *vlp *= (*scale);
      *vup *= (*scale);
    }
  } /* end scaling */
}




/* If matrix scaled, rescale eigenvalues */
static inline 
void invscale_eigvals(val_t *Wstruct, double scale)
{
  int    *size = Wstruct->mp;
  double *vlp  = Wstruct->vlp;
  double *vup  = Wstruct->vup;
  double *W    = Wstruct->W;
  double invscale = 1.0 / scale;
  int    IONE = 1;

  if (scale != 1.0) {  /* FP comparison okay */
    *vlp *= invscale;
    *vup *= invscale;
    dscal_(size, &invscale, W, &IONE);
  }
}





/* Wrapper to call LAPACKs DSTEMR for small matrices */
static
int handle_small_cases(char *jobz, char *range, int *np, double  *D,
		       double *E, double *vl, double *vu, int *il,
		       int *iu, int *tryrac, int *mp, double *W, 
		       double *Z, int *ldz, int *Zsupp)
{
  bool   cntval = (jobz[0]  == 'C' || jobz[0]  == 'c');
  bool   onlyW  = (jobz[0]  == 'N' || jobz[0]  == 'n');
  bool   wantZ  = (jobz[0]  == 'V' || jobz[0]  == 'v');
  int    n          = *np;
  int    lwork, *iwork, liwork, info, MINUSONE=-1;
  double *work;
  double cnt;

  if (onlyW) {
    lwork  = 12*n;
    liwork =  8*n;
  } else if (wantZ || cntval) {
    lwork  = 18*n;
    liwork = 10*n;
  } 

  work = (double *) malloc( lwork  * sizeof(double));
  assert(work != NULL);
  
  iwork = (int *)   malloc( liwork * sizeof(int));
  assert(iwork != NULL);

  if (cntval) {
    dstemr_("V", "V", np, D, E, vl, vu, il, iu, mp, W, &cnt,
	    ldz, &MINUSONE, Zsupp, tryrac, work, &lwork, iwork,
	    &liwork, &info);
    assert(info == 0);
  
    *mp = (int) cnt; 
    free(work); free(iwork);
    return(0);
  }

  dstemr_(jobz, range, np, D, E, vl, vu, il, iu, mp, W, Z,
	  ldz, np, Zsupp, tryrac, work, &lwork, iwork,
	  &liwork, &info);
  assert(info == 0);
  
  free(work);
  free(iwork);

  return(0);
}




/* Refines the eigenvalue to high relative accuracy with
 * respect to the input matrix;
 * Note: In principle this part could be multithreaded too,
 * but it will only rarely be called and not much work
 * is involved */
static 
void refine_to_highrac(double *D, double *E2, in_t *Dstruct, 
		       val_t *Wstruct, tol_t *tolstruct)
{
  int              n      = Dstruct->n;
  int              nsplit = Dstruct->nsplit;
  int    *restrict isplit = Dstruct->isplit;
  double           spdiam = Dstruct->spdiam;
  int              *mp    = Wstruct->mp;
  double *restrict W      = Wstruct->W;
  double *restrict Werr   = Wstruct->Werr;
  int    *restrict Windex = Wstruct->Windex;
  int    *restrict iblock = Wstruct->iblock;
  double           pivmin = tolstruct->pivmin; 
  double           tol    = 4 * DBL_EPSILON; 
  
  double           *work;
  int              *iwork;
  int              i, j;
  int              ilow, iupp, offset, info;
  int              begin,  end,  nbl;
  int              Wbegin, Wend, mbl;

  work  = (double *) malloc( 2*n * sizeof(double) );
  assert (work != NULL);
  iwork = (int *)    malloc( 2*n * sizeof(int)    );
  assert (iwork != NULL);

  begin  = 0;
  Wbegin = 0;
  for (j=0; j<nsplit; j++) {
    
    end = isplit[j] - 1;
    nbl = end - begin + 1;
    mbl = 0;
    
    /* Find eigenvalues in block */
    ilow = n;
    iupp = 1;
    for (i=Wbegin; i<*mp; i++) {
      if (iblock[i] == j+1)  {
	mbl++;
	ilow = imin(ilow, Windex[i]);
	iupp = imax(iupp, Windex[i]);
      } else {
	break;
      }
    }

    /* If no eigenvalues for process in block continue */
    if (mbl == 0) {
      begin  = end  + 1;
      continue;
    }

    Wend = Wbegin + mbl - 1;
    offset  = ilow - 1;

    dlarrj_(&nbl, &D[begin], &E2[begin], &ilow, &iupp, &tol,
	    &offset, &W[Wbegin], &Werr[Wbegin], work, iwork, &pivmin,
	    &spdiam, &info);
    assert(info == 0);
    
    begin  = end  + 1;
    Wbegin = Wend + 1;
  } /* end j */
  
  free(work);
  free(iwork);
}




static inline void sort_eigenpairs(val_t *Wstruct, vec_t *Zstruct,
				   double *work)
{
  int n                   = Wstruct->n;
  int m                   = Wstruct->m;
  double *restrict W      = Wstruct->W;
  int    *restrict Windex = Wstruct->Windex;
  int              ldz    = Zstruct->ldz;
  double *restrict Z      = Zstruct->Z;
  int    *restrict Zsupp  = Zstruct->Zsupp;
  int    *restrict Zindex = Zstruct->Zindex;

  int           i;
  sort_struct_t *array;
  bool          sorted;
  double        tmp;
  int           itmp1, itmp2;

  array = (sort_struct_t *) malloc( m*sizeof(sort_struct_t) );

  for (i=0; i<m; i++) {
    array[i].lambda = W[i]; 
    array[i].ilocal = Windex[i];
    array[i].iblock = 0;
    array[i].ind    = Zindex[i];
  }

  /* Sort according to Zindex */
  qsort(array, m, sizeof(sort_struct_t), cmpb);

  for (i=0; i<m; i++) {
    W[i]      = array[i].lambda; 
    Windex[i] = array[i].ilocal;
  }

  /* Make sure that sorted correctly; ineffective implementation,
   * but usually no or very little swapping should be done here */
  sorted = false;
  while (sorted == false) {
    sorted = true;
    for (i=0; i<m-1; i++) {
      if (W[i] > W[i+1]) {
	sorted = false;
	/* swap eigenvalue */
	tmp    = W[i];
	W[i]   = W[i+1];
	W[i+1] = tmp;
	/* swap eigenvalue support */
	itmp1 = Zsupp[2*i];
	Zsupp[2*i] = Zsupp[2*(i+1)];
	Zsupp[2*(i+1)] = itmp1;
	itmp2 = Zsupp[2*i + 1];
	Zsupp[2*i + 1] = Zsupp[2*(i+1) + 1];
	Zsupp[2*(i+1) + 1] = itmp2;
	/* swap eigenvector (should only copy support) */
	memcpy(work, &Z[i*ldz], n*sizeof(double));
	memcpy(&Z[i*ldz], &Z[(i+1)*ldz], n*sizeof(double));
	memcpy(&Z[(i+1)*ldz], work, n*sizeof(double));
      }
    }
  } /* end while */

  free(array);
}




/* Compare function for using qsort() on an array
 * of doubles */
static int cmpa( const void *a1, const void *a2 )
{
  double arg1 = *(double *)a1;
  double arg2 = *(double *)a2;

  if ( arg1 < arg2 ) {
    return(-1);
  } else {
    return(1);
  }
}



/* Compare function for using qsort() on an array of 
 * sort_structs */
static 
int cmpb( const void *a1, const void *a2 )
{
  sort_struct_t *arg1, *arg2;

  arg1 = (sort_struct_t *) a1;
  arg2 = (sort_struct_t *) a2;

  /* within block local index decides */
  if ( arg1->ind < arg2->ind ) 
    return(-1);
  else
    return(1);
}




/* Fortran function prototype */
void mrrr_(char *jobz, char *range, int *n, double  *D,
	   double *E, double *vl, double *vu, int *il, int *iu,
	   int *tryrac, int *m, double *W, double *Z, int *ldz, 
	   int *Zsupp, int* info)
{
  *info = mrrr(jobz, range, n, D, E, vl, vu, il, iu, tryrac, 
	       m, W, Z, ldz, Zsupp);
}
