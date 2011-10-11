/* Computing all or a subset of eigenvectors, given the eigenvalues
 * of a tridiagonal matrix T preprocessed by 'mrrr_val'.
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
#include <math.h>
#include <assert.h>
#include <pthread.h> 

#include "global.h"
#include "mrrr.h"  
#include "counter.h"
#include "rrr.h"
#include "queue.h"
#include "tasks.h"
#include "structs.h"


static inline workQ_t *create_workQ();
static inline void destroy_workQ(workQ_t*);
static inline aux3_t *create_aux3(int, int, counter_t*, workQ_t*, 
				  in_t*, val_t*, vec_t*, tol_t*);
static void *empty_workQ(void*);
static inline void init_workQ(workQ_t*, in_t*, val_t*);
static inline void init_zindex(in_t*, val_t*, vec_t*);
static int cmp(const void*, const void*);



int mrrr_vec(int nthreads, in_t *Dstruct, val_t *Wstruct, 
	     vec_t *Zstruct, tol_t *tolstruct)
{
  int            n = Dstruct->n;
  int            m = Wstruct->m;
  double         *Wshifted;
  int            i, info;
  pthread_t      *threads;
  pthread_attr_t attr;
  void           *status;
  aux3_t         *arg;
  workQ_t        *workQ;
  counter_t      *num_left;

  Wshifted = (double *) malloc( n * sizeof(double) );
  assert(Wshifted != NULL);

  Wstruct->Wshifted = Wshifted;

  threads = (pthread_t *) malloc(nthreads * sizeof(pthread_t));
  assert(threads != NULL);

  /*  Initialize index vector for eigenvectors */
  init_zindex(Dstruct, Wstruct, Zstruct);

  /* Create work queue, counter, and threads to empty work queue */
  workQ      = create_workQ( );
  num_left   = PMR_create_counter(m);

  threads[0] = pthread_self();
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

  /* Create nthreads-1 additional threads */
  for (i=1; i<nthreads; i++) {

    arg = create_aux3(i, nthreads, num_left, workQ, Dstruct, 
		      Wstruct, Zstruct, tolstruct);

    info = pthread_create(&threads[i], &attr, empty_workQ, 
			  (void *) arg);
    assert(info == 0);
  }
  
  /* Initialize the work queue with tasks */
  init_workQ(workQ, Dstruct, Wstruct);

  arg = create_aux3(0, nthreads, num_left, workQ, Dstruct, 
		    Wstruct, Zstruct, tolstruct);

  status = empty_workQ((void *) arg);
  assert(status == NULL);
 
  /* Join all the worker thread */
  for (i=1; i<nthreads; i++) {
    info = pthread_join(threads[i], &status);
    assert(info == 0 && status == NULL);
  }

  /* Clean up and return */
  free(Wshifted);
  free(threads);
  pthread_attr_destroy(&attr);
  destroy_workQ(workQ);
  PMR_destroy_counter(num_left);
  
  return(0);
}
  



static inline workQ_t *create_workQ()
{
  workQ_t *wq;

  wq = (workQ_t *) malloc( sizeof(workQ_t) );

  wq->r_queue = PMR_create_empty_queue();
  wq->s_queue = PMR_create_empty_queue();
  wq->c_queue = PMR_create_empty_queue();

  return(wq);
}




static inline void destroy_workQ(workQ_t *wq)
{
  PMR_destroy_queue(wq->r_queue);
  PMR_destroy_queue(wq->s_queue);
  PMR_destroy_queue(wq->c_queue);
  free(wq);
}




static inline 
aux3_t *create_aux3(int tid, int nthreads, counter_t *num_left, 
		    workQ_t *workQ, in_t *Dstruct, 
		    val_t *Wstruct, vec_t *Zstruct, 
		    tol_t *tolstruct)
{
  aux3_t *arg;

  arg = (aux3_t *) malloc( sizeof(aux3_t) );
  assert(arg != NULL);

  arg->tid       = tid;
  arg->nthreads  = nthreads;
  arg->num_left  = num_left;
  arg->workQ     = workQ;
  arg->Dstruct   = Dstruct;
  arg->Wstruct   = Wstruct;
  arg->Zstruct   = Zstruct;
  arg->tolstruct = tolstruct;

  return(arg);
}




static inline 
void retrieve_aux3(aux3_t *arg, int *tid, int *nthreads, 
		   counter_t **num_left, workQ_t **workQ, 
		   in_t **Dstruct, val_t **Wstruct, 
		   vec_t **Zstruct, tol_t **tolstruct)
{
  *tid       = arg->tid;
  *nthreads  = arg->nthreads;
  *num_left  = arg->num_left;
  *workQ     = arg->workQ;
  *Dstruct   = arg->Dstruct;
  *Wstruct   = arg->Wstruct;
  *Zstruct   = arg->Zstruct;
  *tolstruct = arg->tolstruct;

  free(arg);
}





/*
 * Processes all the tasks put in the work queue.
 */
static void *empty_workQ(void *argin)
{
  /* input arguments */
  int       tid;
  int       nthreads;
  counter_t *num_left;
  workQ_t   *workQ;
  in_t      *Dstruct;
  val_t     *Wstruct;
  vec_t     *Zstruct;
  tol_t     *tolstruct;
  int       n;

  /* others */
  task_t    *task;
  double    *work;
  int       *iwork;

  /* retrieve necessary arguments from structures */
  retrieve_aux3((aux3_t *) argin, &tid, &nthreads, &num_left,
		&workQ, &Dstruct, &Wstruct, &Zstruct, &tolstruct);
  
  n = Wstruct->n;

  /* max. needed double precision work space: dlar1v */
  work      = (double *) malloc( 4*n * sizeof(double) );
  assert(work != NULL);

  /* max. needed double precision work space: dlarrb */
  iwork     = (int *)    malloc( 2*n * sizeof(int)    );
  assert(iwork != NULL);

  /* While loop to empty the work queue */
  while (PMR_get_counter_value(num_left) > 0) {
    
    task = PMR_remove_task_at_front(workQ->r_queue);
    if (task != NULL) {
      assert(task->flag == REFINEMENT_TASK_FLAG);

      PMR_process_r_task((refine_t *) task->data, tid, Wstruct, 
			 tolstruct, work, iwork);
      free(task);
      continue;
    }
    
    task = PMR_remove_task_at_front(workQ->s_queue);
    if ( task != NULL ) {
      assert(task->flag == SINGLETON_TASK_FLAG);

      PMR_process_s_task((singleton_t *) task->data, tid, num_left, 
			 workQ, Wstruct, Zstruct, tolstruct, work, 
			 iwork);
      free(task);
      continue;
    }
    
    task = PMR_remove_task_at_front(workQ->c_queue);
    if ( task != NULL ) {
      assert(task->flag == CLUSTER_TASK_FLAG);

      PMR_process_c_task((cluster_t *) task->data, tid, nthreads, 
			 num_left, workQ, Wstruct, Zstruct, tolstruct, 
			 work, iwork);
      free(task);
      continue;
    }

  } /* end while */

  free(work);
  free(iwork);

  return(NULL);
}




static inline void init_workQ(workQ_t *workQ, in_t *Dstruct, 
			      val_t *Wstruct)
{
  double *restrict D      = Dstruct->D;
  double *restrict L      = Dstruct->E;
  int    *restrict isplit = Dstruct->isplit;
  int              m      = Wstruct->m;
  double           *vlp   = Wstruct->vlp;
  double *restrict W      = Wstruct->W;
  double *restrict Werr   = Wstruct->Werr;
  int    *restrict iblock = Wstruct->iblock;

  int              j;
  int              begin ,  end;
  int              Wbegin, Wend;
  int              nbl;
  double           sigma;
  rrr_t            *RRR;
  double           lgap;
  task_t           *task;

  /* For every unreducible block of the matrix create a task
   * and put it in the queue */
  begin  = 0;
  Wbegin = 0;

  for (j=0; j<Dstruct->nsplit; j++) {

    end   = isplit[j] - 1;
    sigma = L[end];

    Wend = Wbegin-1;
    while (Wend < m-1 && iblock[Wend + 1] == j+1) {
	Wend++;
    }

    if (Wend < Wbegin) {
      begin = end + 1; 
      continue;
    }

    nbl = end - begin  + 1;

    RRR = PMR_create_rrr(&D[begin], &L[begin], NULL, NULL, nbl, -1);
  
    if (nbl == 1) {
      /* To make sure that RRR is freed when s-task is processed */
      PMR_increment_rrr_dependencies(RRR);
      PMR_set_parent_processed_flag(RRR);

      task = PMR_create_s_task(Wbegin, Wbegin, 1, begin, end, Wbegin, 
			       Wend, 0, 0, RRR);
      PMR_insert_task_at_back(workQ->s_queue, task);
      
      begin  = end  + 1;
      Wbegin = Wend + 1;
      continue;
    }

    lgap = fmax(0.0, (W[Wbegin]+sigma) - Werr[Wbegin] - (*vlp) );

    task = PMR_create_c_task(Wbegin, Wend, 0, begin, end, Wbegin, 
			     Wend, 0, lgap, RRR);
  
    PMR_insert_task_at_back(workQ->c_queue, task);

    begin  = end  + 1;
    Wbegin = Wend + 1;
  }
  /* end of loop over block */
}




static inline 
void init_zindex(in_t *Dstruct, val_t *Wstruct, vec_t *Zstruct)
{
  double *restrict L      = Dstruct->E;
  int    *restrict isplit = Dstruct->isplit;
  int              m      = Wstruct->m;
  double *restrict W      = Wstruct->W;
  int    *restrict Windex = Wstruct->Windex;
  int    *restrict iblock = Wstruct->iblock;
  int    *restrict Zindex = Zstruct->Zindex;

  int           i, j;
  int           iblk, ishift; 
  double        sigma;
  sort_struct_t *array;

  array = (sort_struct_t *) malloc(m*sizeof(sort_struct_t));
  assert(array != NULL);

  for (i=0; i<m; i++) {
    iblk          = iblock[i];
    ishift        = isplit[iblk-1] - 1;
    sigma         = L[ishift];
    array[i].lambda = W[i] + sigma; 
    array[i].ilocal = Windex[i];
    array[i].iblock = iblk;
    array[i].ind    = i;
  }

  qsort(array, m, sizeof(sort_struct_t), cmp);
  
  for (i=0; i<m; i++) {
    j = array[i].ind;
    Zindex[j] = i;
  }

  free(array);
}





/* Compare function for using qsort() on an array of sort_structs */
static 
int cmp(const void *a1, const void *a2)
{
  sort_struct_t *arg1, *arg2;

  arg1 = (sort_struct_t *) a1;
  arg2 = (sort_struct_t *) a2;

  /* Within block local index decides */
  if ( arg1->iblock == arg2->iblock ) {
    return (arg1->ilocal - arg2->ilocal);
  } else {
    if ( arg1->lambda < arg2->lambda ) {
      return(-1);
    } else if ( arg1->lambda > arg2->lambda ) {
      return(1);
    } else {
      if (arg1->ilocal < arg2->ilocal) return(-1);
      else return(1);
    }
  }

}
