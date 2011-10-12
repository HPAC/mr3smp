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

#ifndef TTASK_H
#define TTASK_H

#include <stdlib.h>
#include <semaphore.h>
#include "rrr.h"
#include "queue.h"
#include "counter.h"
#include "structs.h"


#define SINGLETON_TASK_FLAG  0
#define CLUSTER_TASK_FLAG    1
#define REFINEMENT_TASK_FLAG 2


typedef struct {
  int        begin;
  int        end;
  int        depth;
  int        bl_begin;
  int        bl_end;
  int        bl_W_begin;
  int        bl_W_end;
  double     bl_spdiam;
  double     lgap;
  rrr_t      *RRR;
} cluster_t, singleton_t;

typedef struct {
  int        taskcount;
  //----------------------------------------
  //pthread_mutex_t    taskmutex;
  counter_t  *counter;
  //----------------------------------------

  cluster_t  *cl;
  int        nthreads;
  counter_t  *num_left;
  workQ_t    *workQ;
  rrr_t      *RRR;
  vec_t      *Zstruct;
	int 			 num_tasks;
	int				 chunk;
} subtasks_t;

typedef struct {
  int       begin;
  double    *D;
  double    *DLL;
  int       p;
  int       q;
  int       bl_size;
  double    bl_spdiam;
  int       producer;
  subtasks_t *sts;
} refine_t;


task_t *PMR_create_c_task(int cl_begin, int cl_end, int depth, 
			  int bl_begin, int bl_end, int bl_W_begin, 
			  int bl_W_end, double bl_spdiam, double lgap,
			  rrr_t *RRR);

task_t *PMR_create_s_task(int begin, int end, int depth, 
			  int bl_begin, int bl_end, int bl_W_begin, 
			  int bl_W_end, double bl_spdiam, double lgap,
			  rrr_t *RRR);

task_t *PMR_create_r_task(int rf_begin, double *D, double *DLL, int p, 
			  int q, int bl_size, double bl_spdiam, 
                          int producer, subtasks_t *sts);

int PMR_process_c_task(cluster_t *cl, int tid, int nthreads, 
		       counter_t *num_left, workQ_t *workQ, 
		       val_t *Wstruct, vec_t *Zstruct, 
		       tol_t *tolstruct, double *work, int *iwork);

int PMR_process_s_task(singleton_t *sng, int tid, counter_t *num_left, 
		       workQ_t *workQ, val_t *Wstruct, vec_t *Zstruct, 
		       tol_t *tolstruct, double *work, int *iwork);

int PMR_process_r_task(refine_t *rf, int tid, val_t *Wstruct, 
		       tol_t *tolstruct, double *work, int *iwork);

int PMR_process_r_queue(int tid, workQ_t *workQ, val_t *Wstruct, 
			tol_t *tolstruct, double *work, int *iwork);

int PMR_create_subtasks(cluster_t *cl, int tid, int nthreads, 
                    counter_t *num_left, workQ_t *workQ, rrr_t *RRR, 
                    val_t *Wstruct, vec_t *Zstruct, tol_t *tolstruct,
                    double *work, int *iwork);

void PMR_destroy_task(task_t *task);

#endif
