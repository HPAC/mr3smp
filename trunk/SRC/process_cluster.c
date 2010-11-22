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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <semaphore.h>
#include "global.h"
#include "mrrr.h"
#include "counter.h"
#include "queue.h"
#include "tasks.h"
#include "rrr.h"
#include "structs.h"


static inline 
rrr_t* compute_new_rrr(cluster_t *cl, int tid, val_t *Wstruct, 
		       vec_t *Zstruct, tol_t *tolstruct, double *work, 
		       int *iwork);

static inline 
int refine_eigvals(cluster_t *cl, int tid, int nthreads, 
		   counter_t *num_left, workQ_t *workQ, rrr_t *RRR, 
		   val_t *Wstruct, vec_t *Zstruct, tol_t *tolstruct, 
		   double *work, int *iwork);

static inline 
int create_subtasks(cluster_t *cl, int tid, int nthreads, 
		    counter_t *num_left, workQ_t *workQ, rrr_t *RRR, 
		    val_t *Wstruct, vec_t *Zstruct, tol_t *tolstruct,
		    double *work, int *iwork);





/* Processing a cluster */
int PMR_process_c_task(cluster_t *cl, int tid, int nthreads, 
		       counter_t *num_left, workQ_t *workQ, 
		       val_t *Wstruct, vec_t *Zstruct, tol_t *tolstruct, 
		       double *work, int *iwork)
{
  /* Inputs */
  int depth    = cl->depth;
  int m        = Wstruct->m;

  /* Others */
  int   info;
  rrr_t *RRR;

  assert(depth < m);

  RRR = compute_new_rrr(cl, tid, Wstruct, Zstruct, tolstruct, 
			work, iwork);

  info = refine_eigvals(cl, tid, nthreads, num_left, workQ, RRR, 
			Wstruct, Zstruct, tolstruct, work, iwork);
  assert(info == 0);
  
  info = create_subtasks(cl, tid, nthreads, num_left, workQ, RRR, 
			 Wstruct, Zstruct, tolstruct, work, iwork);
  assert(info == 0);

  return(0);
}





static inline 
rrr_t* compute_new_rrr(cluster_t *cl, int tid, val_t *Wstruct, 
		       vec_t *Zstruct, tol_t *tolstruct, 
		       double *work, int *iwork)
{
  int    cl_begin    = cl->begin;
  int    cl_end      = cl->end;
  int    cl_size     = cl_end - cl_begin + 1;
  int    depth       = cl->depth;
  int    bl_begin    = cl->bl_begin;
  int    bl_end      = cl->bl_end;
  int    bl_size     = bl_end - bl_begin + 1;
  double bl_spdiam   = cl->bl_spdiam;
  rrr_t  *RRR_parent = cl->RRR;

  double *restrict Werr     = Wstruct->Werr;
  double *restrict Wgap     = Wstruct->Wgap;
  int    *restrict Windex   = Wstruct->Windex;
  double *restrict Wshifted = Wstruct->Wshifted;
  double           RQtol    = tolstruct->RQtol;
  double           pivmin   = tolstruct->pivmin;

  /* Others */
  double *restrict D,         *restrict L;
  double *restrict DL,        *restrict DLL;
  double *restrict D_parent,  *restrict L_parent;
  double *restrict DL_parent, *restrict DLL_parent;

  int    i, k, p, info;
  double tmp;
  int    offset, IONE=1;
  double lgap, rgap, tau, fudge, savegap;
  rrr_t  *RRR;

  /* Allocate memory for new representation for cluster */
  DL  = (double *) malloc(bl_size * sizeof(double));
  assert(DL != NULL);
  
  DLL = (double *) malloc(bl_size * sizeof(double));
  assert(DLL != NULL);

  if (depth == 0) {

    D = RRR_parent->D;
    L = RRR_parent->L;

  } else {
    
    D = (double *) malloc (bl_size *sizeof(double));
    assert(D != NULL);
    
    L = (double *) malloc (bl_size *sizeof(double));
    assert(L != NULL);

    /* Recompute DL and DLL of parent */
    D_parent = RRR_parent->D;
    L_parent = RRR_parent->L;
    for (i=0; i<bl_size-1; i++) {
	tmp    = D_parent[i]*L_parent[i];
	DL[i]  = tmp;
	DLL[i] = tmp*L_parent[i];
    }
    DL_parent  = DL;
    DLL_parent = DLL;
    
    /* Shift as close as possible refine extremal eigenvalues */
    for (k=0; k<2; k++) {
      if (k == 0) {
	p              = Windex[cl_begin];
	savegap        = Wgap[cl_begin];
	Wgap[cl_begin] = 0.0;
      } else {
	p              = Windex[cl_end];
	savegap        = Wgap[cl_end];
	Wgap[cl_end]   = 0.0;
      }
      
      offset  = Windex[cl_begin] - 1;
      
      dlarrb_(&bl_size, D_parent, DLL_parent, &p, &p, &RQtol,
	      &RQtol, &offset, &Wshifted[cl_begin], &Wgap[cl_begin],
	      &Werr[cl_begin], work, iwork, &pivmin, &bl_spdiam,
	      &bl_size, &info);
      assert(info == 0);
      
      if (k == 0) {
	Wgap[cl_begin] = fmax(0, (Wshifted[cl_begin+1]-Werr[cl_begin+1])
			       - (Wshifted[cl_begin]+Werr[cl_begin]) );
      } else {
	Wgap[cl_end]   = savegap;
      }
    } /* end k */

    /* calculate right and left gaps */
    lgap = cl->lgap;
    rgap = Wgap[cl_end];
    
    /* Compute new RRR and store it in D and L */
    dlarrf_(&bl_size, D_parent, L_parent, DL_parent,
	    &IONE, &cl_size, &Wshifted[cl_begin], &Wgap[cl_begin],
	    &Werr[cl_begin], &bl_spdiam, &lgap, &rgap,
	    &pivmin, &tau, D, L, work, &info);
    assert(info == 0);
    
    /* update shift and store it */
    tmp = L_parent[bl_size-1] + tau;
    L[bl_size-1] = tmp;
    
    /* update shifted eigenvalues and fudge errors */
    for (k=cl_begin; k<=cl_end; k++) {
      fudge  = 3.0 * DBL_EPSILON * fabs( Wshifted[k] );
      Wshifted[k] -= tau;
      fudge += 4.0 * DBL_EPSILON * fabs( Wshifted[k] );
      Werr[k] += fudge;
    } 
  } 
  /* D, L is the RRR of the current depth level now */

  /* Compute corresponding D*L and D*L*L */
  for (i=0; i<bl_size-1; i++) {
    tmp    = D[i]*L[i];
    DL[i]  = tmp;
    DLL[i] = tmp*L[i];
  }
    
  RRR = PMR_reset_rrr(RRR_parent, D, L, DL, DLL, bl_size, depth);

  /* Prevend freeing before cluster task processed */
  PMR_increment_rrr_dependencies(RRR);

  return(RRR);
}



static inline 
int refine_eigvals(cluster_t *cl, int tid, int nthreads, 
		   counter_t *num_left, workQ_t *workQ, rrr_t *RRR, 
		   val_t *Wstruct, vec_t *Zstruct, tol_t *tolstruct, 
		   double *work, int *iwork)
{
  int    cl_begin    = cl->begin;
  int    cl_end      = cl->end;
  int    cl_size     = cl_end - cl_begin + 1;
  int    depth       = cl->depth;
  int    bl_begin    = cl->bl_begin;
  int    bl_end      = cl->bl_end;
  int    bl_size     = bl_end - bl_begin + 1;
  int    bl_W_begin  = cl->bl_W_begin;
  int    bl_W_end    = cl->bl_W_end;
  int    mbl         = bl_W_end - bl_W_begin + 1;
  double bl_spdiam   = cl->bl_spdiam;

  double *restrict D        = RRR->D;
  double *restrict L        = RRR->L;
  double *restrict DLL      = RRR->DLL;
  double *restrict W        = Wstruct->W;
  double *restrict Werr     = Wstruct->Werr;
  double *restrict Wgap     = Wstruct->Wgap;
  int    *restrict Windex   = Wstruct->Windex;
  double *restrict Wshifted = Wstruct->Wshifted;
  double *restrict gersch   = Wstruct->gersch;
  double           rtol1    = tolstruct->rtol1;
  double           rtol2    = tolstruct->rtol2;
  double           pivmin   = tolstruct->pivmin;

  int    i, j, info, offset;
  double gl, gu, sigma, savegap;
  int    nleft, chunk, own_part, others_part;
  int    rf_begin, rf_end, p, q;
  int    num_tasks, count;
  sem_t  sem;
  task_t *task;

  if (depth == 0) {

    gl = gersch[2*bl_begin    ];
    gu = gersch[2*bl_begin + 1];
    for (i = bl_begin+1; i<bl_end; i++) {
      gl = fmin(gl, gersch[2*i    ]);
      gu = fmax(gu, gersch[2*i + 1]);
    }
    cl->bl_spdiam = gu - gl;

    sigma = L[bl_size-1];

    memcpy(&Wshifted[cl_begin], &W[cl_begin], mbl*sizeof(double));
    for (i=0; i<mbl; i++) {
      W[cl_begin + i] += sigma;
    }

  } else {

    /* Determine if refinement should be split into tasks */
    nleft = PMR_get_counter_value(num_left);
    own_part = (int) fmax( ceil( (double) nleft / nthreads ),
			   MIN_BISEC_CHUNK);
    
    if (own_part < cl_size) {
      /* Parallize refinement the refinement by creating tasks */
      num_tasks   = (int) ceil((double) cl_size / own_part) - 1; /*>1*/

      /* Fudge number of tasks to be cerated in case many threads 
       * are involved to get better load balancing; my first version 
       * had gradual smaller task for even better load balancing, 
       * but I thought through the absolute splitting one could get rid 
       * of this slightly ugly fudging */
      if (num_tasks >= (int) ceil((1.0-(1.0/NUM_SOCKETS))*nthreads)) {
	own_part  = (int) fmax( ceil( (double) own_part / 2.0 ),
				MIN_BISEC_CHUNK);
	num_tasks *= 4;
      } 

      others_part = cl_size - own_part;
      chunk       = others_part/num_tasks;
      
      sem_init(&sem, 0, 0);
      rf_begin = cl_begin;
      p        = Windex[cl_begin];
      for (i=0; i<num_tasks; i++) {
	rf_end = rf_begin + chunk - 1;
	q      = p        + chunk - 1;
	
	task = PMR_create_r_task(rf_begin, D, DLL, p, q, bl_size, 
				 bl_spdiam, tid, &sem);

	if (rf_begin <= rf_end)
	  PMR_insert_task_at_back(workQ->r_queue, task);
	else
	  sem_post(&sem); /* case chunk=0 */

	rf_begin = rf_end + 1;
	p        = q      + 1;
      }
      rf_end = cl_end;
      q      = Windex[cl_end];
      offset = Windex[rf_begin] - 1;

      /* Call bisection routine to refine the values */
      if (rf_begin <= rf_end) {
	dlarrb_(&bl_size, D, DLL, &p, &q, &rtol1, &rtol2, &offset, 
		&Wshifted[rf_begin], &Wgap[rf_begin], &Werr[rf_begin],
		work, iwork, &pivmin, &bl_spdiam, &bl_size, &info);
	assert( info == 0 );
      }
      /* Empty all refine tasks before waiting */
      PMR_process_r_queue(tid, workQ, Wstruct, tolstruct, work, 
			  iwork);

      /* Barrier: wait until all created tasks finished */
      count = num_tasks;
      while (count > 0) {
	while (sem_wait(&sem) != 0) { };
	count--;
      }
      sem_destroy(&sem);

      /* Edit right gap at splitting point */
      rf_begin = cl_begin;
      for (i=0; i<num_tasks; i++) {
	rf_end = rf_begin + chunk - 1;
	
	Wgap[rf_end] = Wshifted[rf_end + 1] - Werr[rf_end + 1]
	               - Wshifted[rf_end] - Werr[rf_end];
      
	rf_begin = rf_end + 1;
      }

    } else {
      /* do refinement of cluster without creating tasks */
      
      /* p and q are local (within block) indices of
       * the first/last eigenvalue of the cluster */
      p = Windex[cl_begin];
      q = Windex[cl_end];
      offset = Windex[cl_begin] - 1;
      
      if (p == q) {
	savegap = Wgap[cl_begin];
	Wgap[cl_begin] = 0.0;
      }  
      
      /* call bisection routine to refine the values to rel. acc. */
      dlarrb_(&bl_size, D, DLL, &p, &q, &rtol1, &rtol2, &offset, 
	      &Wshifted[cl_begin], &Wgap[cl_begin], &Werr[cl_begin],
	      work, iwork, &pivmin, &bl_spdiam, &bl_size, &info);
      assert( info == 0 );
      
      if (p == q) {
	Wgap[cl_begin] = savegap;
      }  
    } /* end refine with or without creating tasks */
  
    sigma = L[bl_size-1];
    
    /* refined eigenvalues with all shifts applied in W */
    for ( j=cl_begin; j<=cl_end; j++ ) {
      W[j] = Wshifted[j] + sigma;
    }

  } /* end refining eigenvalues */

  return(0);
}




static inline 
int create_subtasks(cluster_t *cl, int tid, int nthreads, 
		    counter_t *num_left, workQ_t *workQ, rrr_t *RRR, 
		    val_t *Wstruct, vec_t *Zstruct, tol_t *tolstruct,
		    double *work, int *iwork)
{
  /* Inputs */
  int    cl_begin    = cl->begin;
  int    cl_end      = cl->end;
  int    cl_size     = cl_end - cl_begin + 1;
  int    depth       = cl->depth;
  int    bl_begin    = cl->bl_begin;
  int    bl_end      = cl->bl_end;
  int    bl_size     = bl_end - bl_begin + 1;
  int    bl_W_begin  = cl->bl_W_begin;
  int    bl_W_end    = cl->bl_W_end;
  double bl_spdiam   = cl->bl_spdiam;

  double *restrict D        = RRR->D;
  double *restrict L        = RRR->L;
  double *restrict Wgap     = Wstruct->Wgap;
  double *restrict Wshifted = Wstruct->Wshifted;
  int              ldz      = Zstruct->ldz;
  double *restrict Z        = Zstruct->Z;
  int    *restrict Zindex   = Zstruct->Zindex;
  
  /* Others */
  int       i, k;
  int       new_first, new_last;
  int       new_size;
  int       new_ftt0, new_ftt1;
  int       cl_first, cl_last;
  double    lgap;
  rrr_t     *RRR_parent;
  task_t    *task;
  int       sn_first=-1, sn_last=-1, sn_size=-1;
  int       max_size;
  bool      task_inserted;
  double    avggap, avggap_factor;

  max_size = fmax(1, PMR_get_counter_value(num_left) /
		     (fmin(depth+1,4)*nthreads) );
  task_inserted = true;
  avggap_factor = (nthreads > 24) ? 1.0 : 2.0;

  new_first = cl_begin;
  for (i=cl_begin; i<=cl_end; i++) {    

    if ( i == cl_end )
      new_last = i;
    else if ( Wgap[i] >= MIN_RELGAP*fabs(Wshifted[i]) )
      new_last = i;
    else
      continue;

    new_size = new_last - new_first + 1;
    
    if (new_size == 1) {
      /* singleton was found */
      
      if (new_first==cl_begin || task_inserted==true) {
	/* initialize new singleton task */
	sn_first = new_first;
	sn_last  = new_first;
	sn_size  = 1;
      } else {
	/* extend singleton task by one */
	sn_last++;
	sn_size++;
      }
      
      /* insert task if ... */
      if (i==cl_end || sn_size>=max_size ||
	  Wgap[i+1] < MIN_RELGAP*fabs(Wshifted[i+1])) {

	if (sn_first == cl_begin) {
	  lgap = cl->lgap;
	} else {
	  lgap = Wgap[sn_first-1];
	}
	
	PMR_increment_rrr_dependencies(RRR);
	
	task = PMR_create_s_task(sn_first, sn_last, depth+1, bl_begin,
				 bl_end, bl_W_begin, bl_W_end, 
				 bl_spdiam, lgap, RRR);


	if (depth > 0 && cl_size < STASK_NOENQUEUE) {
	  PMR_process_s_task((singleton_t *) task->data, tid, num_left, 
			     workQ, Wstruct, Zstruct, tolstruct, work, 
			     iwork);
	  free(task);
	} else {
	  PMR_insert_task_at_back(workQ->s_queue, task);
	}
	  
	task_inserted = true;
      } else {
	task_inserted = false;
      }

    } else {
      /* cluster was found */
      
      if (depth == 0 && new_size > 3) {
	/* split cluster to smaller clusters by an absolute
         * split criterion */
	
	cl_first = new_first;
	cl_last  = new_last;
	avggap   = bl_spdiam / (bl_size-1);

	for (k=new_first+1; k<new_last; k++) {
	  
	  if (k == new_last-1)
	    cl_last = new_last;
	  else if (Wgap[k] > avggap_factor*avggap)
	    cl_last = k;
	  else
	    continue;

	  /* find gap to the left */
	  if (cl_first == cl_begin) {
	    lgap = cl->lgap;
	  } else {
	    lgap = Wgap[cl_first - 1];
	  }

	  RRR_parent = PMR_create_rrr(D, L, NULL, NULL, bl_size, depth);

	  /* Create the task for the cluster and put it in the queue */ 
	  task = PMR_create_c_task(cl_first, cl_last, depth+1, 
				   bl_begin, bl_end, bl_W_begin, 
				   bl_W_end, bl_spdiam, lgap, 
				   RRR_parent);
	  
	  PMR_insert_task_at_back(workQ->c_queue, task);
	  
	  cl_first = k + 1;
	} /* end k */

      } else {
	
	/* find gap to the left */
	if (new_first == cl_begin) {
	  lgap = cl->lgap;
	} else {
	  lgap = Wgap[new_first - 1];
	}
	
	new_ftt0 = Zindex[new_first    ];
	new_ftt1 = Zindex[new_first + 1];
	
	if (depth == 0) {
	  
	  RRR_parent = PMR_create_rrr(D, L, NULL, NULL, bl_size, depth);
	  
	} else {
	  
	  memcpy(&Z[new_ftt0*ldz+bl_begin], D, bl_size * sizeof(double));
	  memcpy(&Z[new_ftt1*ldz+bl_begin], L, bl_size * sizeof(double));
	  
	  RRR_parent = PMR_create_rrr(&Z[new_ftt0*ldz + bl_begin],
				      &Z[new_ftt1*ldz + bl_begin],
				      NULL, NULL, bl_size, depth);
	}
	
	/* Create the task for the cluster and put it in the queue */ 
	task = PMR_create_c_task(new_first, new_last, depth+1, 
				 bl_begin, bl_end, bl_W_begin, 
				 bl_W_end, bl_spdiam, lgap, 
				 RRR_parent);
	
	PMR_insert_task_at_back(workQ->c_queue, task);
      }	
      
      task_inserted = true;
      
    } /* if singleton or cluster found */

    new_first = i + 1;
  }  /* end i */

  /* Set flag in RRR that last singleton is created */
  PMR_set_parent_processed_flag(RRR);
  
  /* Clean up */
  PMR_try_destroy_rrr(RRR);
  free(cl);

  return(0);
} /* end create_subtasks */
