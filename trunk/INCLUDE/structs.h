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

#ifndef SSTRUCTS_H
#define SSTRUCTS_H

#include "global.h"
#include "counter.h"
#include "queue.h"


typedef struct {
  int              n;
  double *restrict D;
  double *restrict E;
  double *restrict E2;
  int              nsplit;
  int    *restrict isplit ;
  double           spdiam;
} in_t;


typedef struct {
  int              n;
  int              m;
  double           *vlp;
  double           *vup;
  int              *ilp;
  int              *iup;
  int              *mp;
  double *restrict W;
  double *restrict Werr;
  double *restrict Wgap;
  int    *restrict Windex;
  int    *restrict iblock;
  double *restrict Wshifted;
  double *restrict gersch;
} val_t;


typedef struct {
  int              ldz;
  int              nz;
  double *restrict Z;
  int    *restrict Zsupp;
  int    *restrict Zindex;
} vec_t;


typedef struct {
  double split;
  double rtol1;
  double rtol2;
  double RQtol;
  double pivmin;
  double bsrtol;
} tol_t;


typedef struct {
  queue_t *r_queue;
  queue_t *s_queue;
  queue_t *c_queue;
} workQ_t;


typedef struct {
  double lambda;
  int    ilocal;
  int    iblock;
  int    ind;
} sort_struct_t;


typedef struct {
  int       tid;
  int       il;
  int       iu;
  int       my_il;
  int       my_iu;
  int       n;
  double    *restrict D;
  double    *restrict E;
  double    *restrict E2;
  double    *restrict W;
  double    *restrict Werr;
  int           *restrict Windex;
  double    *restrict gersch;
  tol_t        *tolstruct;
} aux1_t;


typedef struct {
  int    tid;
  int    nbl;
  int    rf_begin;
  int    rf_end;
  double spdiam;
  double *D;
  double *DE2;
  val_t  *Wstruct;
  tol_t  *tolstruct;
} aux2_t;


typedef struct {
  int tid;
  int nthreads;
  counter_t *num_left;
  workQ_t   *workQ;
  in_t      *Dstruct;
  val_t     *Wstruct;
  vec_t     *Zstruct;
  tol_t     *tolstruct;
} aux3_t;


#endif
