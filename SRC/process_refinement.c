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
 * September 2010, modified by Elmar Peise, September 2011
 * Version 1.1
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
#include <assert.h>
#include <semaphore.h>
#include "global.h"
#include "mrrr.h"
#include "tasks.h"


int PMR_process_r_task(refine_t *rf, int tid, val_t *Wstruct, 
		       tol_t *tolstruct, double *work, int *iwork)
{
  int    rf_begin  = rf->begin;
  double *D        = rf->D;
  double *DLL      = rf->DLL;
  int    p         = rf->p;
  int    q         = rf->q;
  int    bl_size   = rf->bl_size;
  double bl_spdiam = rf->bl_spdiam;
  subtasks_t  *sts = rf->sts;

  double *Wshifted = Wstruct->Wshifted;
  double *Werr     = Wstruct->Werr;
  double *Wgap     = Wstruct->Wgap;
  int    *Windex   = Wstruct->Windex;
  double rtol1     = tolstruct->rtol1;
  double rtol2     = tolstruct->rtol2;
  double pivmin    = tolstruct->pivmin;

  /* Others */
  int info, offset, taskcount, rf_end, i;
  double sigma;
  double *restrict L;
  double *restrict W;

  offset = Windex[rf_begin] - 1;

  /* Bisection to refine the eigenvalues */
  dlarrb_(&bl_size, D, DLL, &p, &q, &rtol1, &rtol2,
	  &offset, &Wshifted[rf_begin], &Wgap[rf_begin], &Werr[rf_begin],
	  work, iwork, &pivmin, &bl_spdiam, &bl_size, &info);
  assert(info == 0);

  taskcount = PMR_decrement_counter(sts->counter, 1);
  
  if (taskcount == 0) {
    L = sts->RRR->L;
    W = Wstruct->W;
    rf_begin = sts->cl->begin;
    for (i=0; i<sts->num_tasks; i++) {
      rf_end = rf_begin + sts->chunk - 1;
      
      Wgap[rf_end] = Wshifted[rf_end + 1] - Werr[rf_end + 1]
	- Wshifted[rf_end] - Werr[rf_end];
      
      rf_begin = rf_end + 1;
    }
    sigma = L[bl_size-1];
    
    /* refined eigenvalues with all shifts applied in W */
    for ( i=sts->cl->begin; i<=sts->cl->end; i++ ) {
      W[i] = Wshifted[i] + sigma;
    }
    
    /* create subtasks */
    info = PMR_create_subtasks(sts->cl, tid, sts->nthreads, sts->num_left, 
			       sts->workQ, sts->RRR, Wstruct, 
			       sts->Zstruct, tolstruct, work, iwork);
    assert(info == 0);

    PMR_destroy_counter(sts->counter); 
    free(sts);
  }
  
  free(rf);

  return(0);
}




int PMR_process_r_queue(int tid, workQ_t *workQ, val_t *Wstruct, 
			tol_t *tolstruct, double *work, int *iwork)
{
  task_t *task;
  int    info;

  while ((task = PMR_remove_task_at_front(workQ->r_queue)) != NULL) {
    
    assert(task->flag == REFINEMENT_TASK_FLAG);
    
    info = PMR_process_r_task((refine_t *) task->data, tid, Wstruct,
			      tolstruct, work, iwork);
    assert(info == 0);

    free(task);
  }
  
  return(0);
}
