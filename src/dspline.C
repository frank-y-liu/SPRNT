/*********************************************************************
#
#  Copyright (C) 2011, 2014 International Business Machines
#
#  Author:  Frank Liu, IBM
#
#  All Rights Reserved. This program and the accompanying materials
#  are made available under the terms of the Eclipse Public License v1.0
#  which accompanies this distribution, and is available at
#  http://www.eclipse.org/legal/epl-v10.html
#
*********************************************************************/

/************************************************************************
  implementation of dspline
************************************************************************/

#include <assert.h>
#include <stdio.h>

#include "dspline.h"

void DSpline::Init(int npt, int k) {
  assert(npt > 0);
  assert(k >= 1 && k <= 5);

  _k = k; // although cubic is preferred
  _nx = npt;

  // allocate memory
  _t.size(_nx + _k);
  _work.size(2 * _nx);
  _q.size((2 * _k - 1) * _nx);
  _bcoef.size(_nx);
}

int DSpline::Construct(double *xx, double *yy) {

  FORTRAN(dbknot)(&xx[0], &_nx, &_k, &_t[0]);
  FORTRAN(dbintk)
  (&xx[0], &yy[0], &_t[0], &_nx, &_k, &_bcoef[0], &_q[0], &_work[0]);

  _xb = xx[0];
  _xe = xx[_nx - 1];

  return 0;
}

double DSpline::EvalFun(double x) {
  const int idx = 0;
  double y;

  assert(x >= _xb && x <= _xe);
  y = FORTRAN(dbvalu)(&_t[0], &_bcoef[0], &_nx, &_k, &idx, &x, &_inbv,
                      &_work[0]);
  return y;
}

double DSpline::EvalDer(double x) {
  const int idx = 1;
  double dy;

  assert(x >= _xb && x <= _xe);
  dy = FORTRAN(dbvalu)(&_t[0], &_bcoef[0], &_nx, &_k, &idx, &x, &_inbv,
                       &_work[0]);
  return dy;
}

int DSpline::EvalFunAndDer(double x, double &y, double &dy) {
  int idx = 0;

  assert(x >= _xb && x <= _xe);

  idx = 0;
  y = FORTRAN(dbvalu)(&_t[0], &_bcoef[0], &_nx, &_k, &idx, &x, &_inbv,
                      &_work[0]);

  idx = 1;
  dy = FORTRAN(dbvalu)(&_t[0], &_bcoef[0], &_nx, &_k, &idx, &x, &_inbv,
                       &_work[0]);

  return 0;
}

// end
