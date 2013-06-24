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
 * September 2010, Version 1.2
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
#include <assert.h>
#include <semaphore.h>
#include "global.h"
#include "queue.h"
#include "counter.h"
#include "tasks.h"



task_t *PMR_create_c_task(int cl_begin, int cl_end, int depth, 
			  int bl_begin, int bl_end, int bl_W_begin, 
			  int bl_W_end, double bl_spdiam, double lgap,
			  rrr_t *RRR)
{
  cluster_t *c;
  task_t    *t;

  c = (cluster_t *) malloc(sizeof(cluster_t));
  assert(c != NULL);
  t = (task_t *) malloc(sizeof(task_t));
  assert(t != NULL);
  
  c->begin      = cl_begin;
  c->end        = cl_end;
  c->depth      = depth;  
  c->bl_begin   = bl_begin;
  c->bl_end     = bl_end;
  c->bl_W_begin = bl_W_begin;
  c->bl_W_end   = bl_W_end;
  c->bl_spdiam  = bl_spdiam;
  c->lgap       = lgap;
  c->RRR        = RRR;

  t->data = (void *) c;
  t->flag = CLUSTER_TASK_FLAG;
  t->next = NULL;
  t->prev = NULL;

  return(t);
}



task_t *PMR_create_s_task(int begin, int end, int depth, int bl_begin, 
			  int bl_end, int bl_W_begin, int bl_W_end, 
			  double bl_spdiam, double lgap, rrr_t *RRR)
{
  singleton_t *s;
  task_t      *t;

  s = (singleton_t *) malloc(sizeof(singleton_t));
  assert(s != NULL);
  t = (task_t *) malloc(sizeof(task_t));
  assert(t != NULL);
  
  s->begin      = begin;
  s->end        = end;
  s->depth      = depth;  
  s->bl_begin   = bl_begin;
  s->bl_end     = bl_end;
  s->bl_W_begin = bl_W_begin;
  s->bl_W_end   = bl_W_end;
  s->bl_spdiam  = bl_spdiam;
  s->lgap       = lgap;
  s->RRR        = RRR;

  t->data = (void *) s;
  t->flag = SINGLETON_TASK_FLAG;
  t->next = NULL;
  t->prev = NULL;

  return(t);
}



task_t *PMR_create_r_task(int rf_begin, double *D, double *DLL, int p, 
			  int q, int bl_size, double bl_spdiam, 
			  int producer, subtasks_t *sts)
{
  refine_t *r;
  task_t   *t;

  r = (refine_t *) malloc(sizeof(refine_t));
  assert(r != NULL);
  t = (task_t *) malloc(sizeof(task_t));
  assert(t != NULL);

  r->begin     = rf_begin;
  r->D         = D;
  r->DLL       = DLL;
  r->p         = p;
  r->q         = q;
  r->bl_size   = bl_size;
  r->bl_spdiam = bl_spdiam;
  r->producer  = producer;
  r->sts       = sts;
  
  t->data = (void *) r;
  t->flag = REFINEMENT_TASK_FLAG;
  t->next = NULL;
  t->prev = NULL;

  return(t);
}



void PMR_destroy_task(task_t *task)
{
  free(task);
}








