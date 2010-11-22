/* Copyright (c) 2010, RWTH Aachen University
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
 * September 2010, Version 1.0
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

#ifndef EEIGTRI_H
#define EEIGTRI_H

#include "global.h"

/* Parallel computation of all or a subset of eigenvalues and 
 * optionally eigenvectors of a symmetric tridiagonal matrix based on 
 * the algorithm of Multiple Relatively Robust Representations (MRRR). 
 * The routine targets multi-core architectures. 
 * The implementation is based on LAPACK's routine 'dstemr'.
 *
 * Function prototype: */

int mrrr(char *jobz, char *range, int *n, double *restrict D,
	 double *restrict E, double *vl, double *vu, int *il, int *iu,
	 bool *tryrac, int *m, double *W, double *Z, int *ldz,
	 int *Zsupp);

/* Arguments:
 * ----------
 *
 * INPUTS: 
 * -------
 * jobz              "N" or "n" - compute only eigenvalues
 *                   "V" or "v" - compute also eigenvectors
 *                   "C" or "c" - count the maximal number of 
 *                                locally computed eigenvectors
 * range             "A" or "a" - all
 *                   "V" or "v" - by interval: (VL,VU]
 *                   "I" or "i" - by index:     IL-IU
 * n                 Matrix size
 * ldz               Leading dimension of eigenvector matrix Z; 
 *                   often equal to matrix size n
 *
 * INPUT + OUTPUT: 
 * ---------------
 * D (double[n])     Diagonal elements of tridiagonal T.
 *                   (On output the array will be overwritten).
 * E (double[n])     Off-diagonal elements of tridiagonal T.
 *                   First n-1 elements contain off-diagonals,
 *                   the last element can have an abitrary value. 
 *                   (On output the array will be overwritten.)
 * vl                If range="V", lower bound of interval
 *                   (vl,vu], on output refined.
 *                   If range="A" or "I" not referenced as input.
 *                   On output the interval (vl,vu] contains ALL
 *                   the computed eigenvalues.
 * vu                If range="V", upper bound of interval
 *                   (vl,vu], on output refined.
 *                   If range="A" or "I" not referenced as input.
 *                   On output the interval (vl,vu] contains ALL
 *                   the computed eigenvalues.
 * il                If range="I", lower index (1-based indexing) of 
 *                   the subset 'il' to 'iu'.
 *                   If range="A" or "V" not referenced as input.
 *                   On output the eigenvalues with index il to iu are 
 *                   computed.
 * iu                If range="I", upper index (1-based indexing) of 
 *                   the subset 'il' to 'iu'.
 *                   If range="A" or "V" not referenced as input.
 *                   On output the eigenvalues with index il to iu are 
 *                   computed.
 * tryrac            0 - do not try to achieve high relative accuracy.
 *                   1 - relative accuracy will be attempted; 
 *                       on output it is set to zero if high relative 
 *                       accuracy is not achieved.
 *
 * OUTPUT: 
 * -------
 * m                 Number of eigenvalues and eigenvectors computed. 
 *                   If jobz="C", 'm' will be set to the number of 
 *                   eigenvalues/eigenvectors that will be computed. 
 * W (double[n])     Eigenvalues
 *                   The first 'm' entries contain the eigenvalues.
 * Z                 Eigenvectors.
 * (double[n*m])     Enough space must be provided to store the
 *                   vectors. 'm' should be bigger or equal 
 *                   to 'n' for range="A" and 'iu-il+1' for range="I".
 *                   For range="V" the minimum 'm' can be queried using
 *                   jobz="C".
 * Zsupp             Support of eigenvectors, which is given by
 * (double[2*n])     Z[2*i] to Z[2*i+1] for the i-th eigenvector
 *                   stored locally (1-based indexing).
 *
 * RETURN VALUE: 
 * -------------
 *                 0 - Success  
 *                 1 - Wrong input parameter
 *                 2 - Misc errors  
 *
 * The Fortran interface takes an additinal integer argument INFO
 * to retrieve the return value. 
 *
 * CALL MRRR(JOBZ, RANGE, N, D, E, VL, VU, IL, IU, TRYRAC, 
 *            M, W, Z, LDZ, ZSUPP, INFO)
 * 
 * CHARACTER        JOBZ, RANGE
 * INTEGER          N, IL, IU, TRYRAC, M, LDZ, ZSUPP(*), INFO
 * DOUBLE PRECISION D(*), E(*), VL, VU, W(*), Z(*,*)
 *
 */


/* Set the number of threads in case OMP_NUM_THREADS is not 
 * specified */
#define DEFAULT_NUM_THREADS    4

/* Set the minumum matrix size for which the multi-threaded code
 * is called. If the matrix is smaller DSTEMR will be called */  
#define DSTEMR_IF_SMALLER     500

/* Set flag to force the use of bisection for the initial eigenvalues 
 * approximation; defualt: false */ 
#define FORCE_BISECTION      false

/* Set flag to switch to bisection if enough parallelism is available,
 * that is if '#cores' > SWITCH_TO_BISEC * '#eigenvalues'/ n;
 * the optimal values is matrix dependent at lies between 8 and 16, 
 * the default is therefore 12 */
#define SWITCH_TO_BISEC       12

/* In case of an SMP system, number of sockets; used to split the 
 * refinement of eigenvalues in smaller tasks if more than or equals
 * (1 - 1/NUM_SOCKETS)*'#threads' threads would be involved anyway;
 * default: 1 */
#define NUM_SOCKETS            1

/* Set flag if Rayleigh Quotient Correction should be used, 
 * which is usually faster; default: true */
#define TRY_RQC              true  

/* Set the maximal allowed element growth for being accepted as 
 * an RRR, that is if max. pivot < MAX_GROWTH * 'spectral diameter'
 * the RRR is accepted; default: 64.0 */
#define MAX_GROWTH           64.0

/* Set the min. relative gap for an eigenvalue to be considered 
 * well separated, that is a singleton; this is a very important 
 * parameter of the computation; default: 10e-3 */
#define MIN_RELGAP           1e-3

/* Set how many iterations should be executed to find the root 
 * representation; default: 6 */
#define MAX_TRY_ROOT          6 

/* Set up the minimum size for the refinement of eigenvalues to 
 * be possibly done in parallel (> 1); default: (very low) 10 */
#define MIN_BISEC_CHUNK      10

/* For any cluster, do not enqueue an S-task if the cluster is 
 * smaller than STASK_NOENQUEUE */
#define STASK_NOENQUEUE      10


/* LAPACK function prototypes of routines that need to be available;
 * Note: type specifier 'extern' does not matter in declaration
 * so here used to mark routines from LAPACK and BLAS libraries */
extern double dlanst_(char*, int*, double*, double*);
extern void   dlarrr_(int*, double*, double*, int*);
extern void   dlarra_(int*, double*, double*, double*, double*, double*,
		      int*, int*, int*);
extern void   dlarrc_(char*, int*, double*, double*, double*, double*,
		      double*, int*, int*, int*, int*);
extern void   dlarrd_(char*, char*, int*, double*, double*, int*, int*,
		      double*, double*, double*, double*, double*, double*,
		      int*, int*, int*, double*, double*, double*, double*,
		      int*, int*, double*, int*, int*);
extern void   dlarrb_(int*, double*, double*, int*, int*, double*,
		      double*, int*, double*, double*, double*, double*,
		      int*, double*, double*, int*, int*);
extern void   dlarnv_(int*, int*, int*, double*);
extern void   dlarrk_(int*, int*, double*, double*, double*, double*,
		      double*, double*, double*, double*, int*);
extern void   dlasq2_(int*, double*, int*);
extern void   dlaebz_(int*, int*, int*, int*, int*, int*,
		      double*, double*, double*, double*, double*, double*,
		      int*, double*, double*, int*, int*, double*,
		      int*, int*);
extern void   dlarrf_(int*, double*, double*, double*, int*, int*, double*,
		      double*, double*, double*, double*, double*, double*,
		      double*, double*, double*, double*, int*);
extern void   dlar1v_(int*, int*, int*, double*, double*, double*, double*,
		      double*, double*, double*, double*, bool*, int*,
		      double*, double*, int*, int*, double*, double*,
		      double*, double*);
extern void   dlarrj_(int*, double*, double*, int*, int*, double*, 
		      int*, double*, double*, double*, int*, double*, 
		      double*, int*);
extern void   dstemr_(char*, char*, int*, double*, double*, double*, 
		      double*, int*, int*, int*, double*, double*, 
		      int*, int*, int*, int*, double*, int*, int*, 
		      int*, int*);

/* BLAS function prototypes - renaming for preventing multithreaded
 * version to be called when linked to a BLAS implementation */
extern void   dscal_(int*, double*, double*, int*);
       void  odscal_(int*, double*, double*, int*);


/* LAPACK routines for wrapper functions of the symmetric and the Hermitian 
 * eigenproblem */ 
extern void dsytrd_(char*, int*, double*, int*, double*, double*, 
		    double*, double*, int*, int*);
extern void dormtr_(char*, char*, char*, int*, int*, double*, int*, 
		    double*, double*, int*, double*, int*, int*);

#ifdef COMPLEX_SUPPORTED
#include <complex.h>
extern void zhetrd_(char*, int*, double complex*, int*, double*, 
		    double*, double complex*, double complex*, int*, 
		    int*);
extern void zunmtr_(char*, char*, char*, int*, int*, double complex*, 
		    int*, double complex*, double complex*, int*, 
		    double complex*, int*, int*);
#endif


/* Wrapper routines for the symmetric and Hermitian eigenproblem
 * Function prototype for the symetric case: */

int dsyeig(char *jobz, char *range, char *uplo, int *n, double *A, 
	   int *lda, double *vl, double *vu, int *il, int *iu, 
	   int *m, double *W, double *Z, int *ldz);

/* Arguments:
 * ----------
 *
 * INPUTS: 
 * -------
 * jobz              "N" - compute only eigenvalues
 *                   "V" - compute also eigenvectors
 * range             "A" - all
 *                   "V" - by interval: (VL,VU]
 *                   "I" - by index:     IL-IU
 * uplo              "L" - Upper triangle of A stored
 *                   "U" - Lower triangle of A stored
 * n                 Order of the matrix A
 * lda               Leading dimension of matrix A; 
 *                   often equal to matrix size n
 * ldz               Leading dimension of eigenvector matrix Z; 
 *                   often equal to matrix size n
 *
 * INPUT + OUTPUT: 
 * ---------------
 * A                 On entry symmetric input matrix, stored 
 * (double[lda*n])   in column major ordering. Depending on the 
 *                   value of 'uplo' only the upper or lower 
 *                   triangular part is referenced.
 *                   (On output the array will be overwritten).
 * vl                If range="V", lower bound of interval
 *                   (vl,vu], on output refined.
 *                   If range="A" or "I" not referenced as input.
 *                   On output the interval (vl,vu] contains ALL
 *                   the computed eigenvalues.
 * vu                If range="V", upper bound of interval
 *                   (vl,vu], on output refined.
 *                   If range="A" or "I" not referenced as input.
 *                   On output the interval (vl,vu] contains ALL
 *                   the computed eigenvalues.
 * il                If range="I", lower index (1-based indexing) of 
 *                   the subset 'il' to 'iu'.
 *                   If range="A" or "V" not referenced as input.
 *                   On output the eigenvalues with index il to iu are 
 *                   computed.
 * iu                If range="I", upper index (1-based indexing) of 
 *                   the subset 'il' to 'iu'.
 *                   If range="A" or "V" not referenced as input.
 *                   On output the eigenvalues with index il to iu are 
 *                   computed.
 *
 * OUTPUT: 
 * -------
 * m                 Number of eigenvalues and eigenvectors computed. 
 * W (double[n])     Eigenvalues
 *                   The first 'm' entries contain the eigenvalues.
 * Z                 Eigenvectors.
 * (double[n*m])     Enough space must be provided to store the
 *                   vectors. 'm' should be bigger or equal 
 *                   to 'n' for range="A" or "V" and 'iu-il+1' for 
 *                   range="I".
 *
 * RETURN VALUE: 
 * -------------
 *                 0 - Success  
 *                 1 - Wrong input parameter
 *                 2 - Misc errors  
 *
 * The Fortran interface takes an additinal integer argument INFO
 * to retrieve the return value. 
 *
 * CALL DSYEIG(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, 
 *             M, W, Z, LDZ, INFO)
 * 
 * CHARACTER        JOBZ, RANGE, UPLO
 * INTEGER          N, IL, IU, M, LDA, LDZ, INFO
 * DOUBLE PRECISION A(*,*), VL, VU, W(*), Z(*,*)
 *
 */
 


/* Function prototype for the symetric case: */

#ifdef COMPLEX_SUPPORTED
int zheeig(char *jobz, char *range, char *uplo, int *np, 
	   double complex *A, int *ldap, double *vlp, double *vup, 
	   int *ilp, int *iup, int *mp, double *W, double complex *Z, 
	   int *ldzp);

/* Arguments:
 * ----------
 *
 * INPUTS: 
 * -------
 * jobz              "N" - compute only eigenvalues
 *                   "V" - compute also eigenvectors
 * range             "A" - all
 *                   "V" - by interval: (VL,VU]
 *                   "I" - by index:     IL-IU
 * uplo              "L" - Upper triangle of A stored
 *                   "U" - Lower triangle of A stored
 * n                 Order of the matrix A
 * lda               Leading dimension of matrix A; 
 *                   often equal to matrix size n
 * ldz               Leading dimension of eigenvector matrix Z; 
 *                   often equal to matrix size n
 *
 * INPUT + OUTPUT: 
 * ---------------
 * A                 On entry Hermitian input matrix, stored 
 * (double complex   in column major ordering. Depending on the
 *  [lda*n])         value of 'uplo' only the upper or lower 
 *                   triangular part is referenced.
 *                   (On output the array will be overwritten).
 * vl                If range="V", lower bound of interval
 *                   (vl,vu], on output refined.
 *                   If range="A" or "I" not referenced as input.
 *                   On output the interval (vl,vu] contains ALL
 *                   the computed eigenvalues.
 * vu                If range="V", upper bound of interval
 *                   (vl,vu], on output refined.
 *                   If range="A" or "I" not referenced as input.
 *                   On output the interval (vl,vu] contains ALL
 *                   the computed eigenvalues.
 * il                If range="I", lower index (1-based indexing) of 
 *                   the subset 'il' to 'iu'.
 *                   If range="A" or "V" not referenced as input.
 *                   On output the eigenvalues with index il to iu are 
 *                   computed.
 * iu                If range="I", upper index (1-based indexing) of 
 *                   the subset 'il' to 'iu'.
 *                   If range="A" or "V" not referenced as input.
 *                   On output the eigenvalues with index il to iu are 
 *                   computed.
 *
 * OUTPUT: 
 * -------
 * m                 Number of eigenvalues and eigenvectors computed. 
 * W (double[n])     Eigenvalues
 *                   The first 'm' entries contain the eigenvalues.
 * Z                 Eigenvectors.
 * (double complex   Enough space must be provided to store the
 *  [n*m])           vectors. 'm' should be bigger or equal 
 *                   to 'n' for range="A" or "V" and 'iu-il+1' for 
 *                   range="I".
 *
 * RETURN VALUE: 
 * -------------
 *                 0 - Success  
 *                 1 - Wrong input parameter
 *                 2 - Misc errors  
 *
 * The Fortran interface takes an additinal integer argument INFO
 * to retrieve the return value. 
 *
 * CALL ZHEEIG(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, 
 *             M, W, Z, LDZ, INFO)
 * 
 * CHARACTER        JOBZ, RANGE, UPLO
 * INTEGER          N, IL, IU, M, LDA, LDZ, INFO
 * DOUBLE PRECISION VL, VU, W(*)
 * COMPLEX*16       A(*,*), Z(*,*)
 *
 */

#endif

#endif







