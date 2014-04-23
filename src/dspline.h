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

/********************************************************
  interface code to dbsplin package
********************************************************/

#ifndef _DSPLINE_H
#define _DSPLINE_H

#include "flexvec.h"

#ifdef LINUX
#define FORTRAN(a) a##_
#else
#define FORTRAN(a) a
#endif

// declare external library
extern "C" void FORTRAN(dbknot)(...);
extern "C" void FORTRAN(dbintk)(...);
extern "C" double FORTRAN(dbvalu)(...);

const int SM_DFLT_FIX = 100;

class DSpline {
 protected:
  int         _nx;
  int         _k;
  
  double      _xb;
  double      _xe;

  int         _inbv;

  FlexVec<double, SM_DFLT_FIX>   _t;
  FlexVec<double, SM_DFLT_FIX>   _bcoef;
  FlexVec<double, SM_DFLT_FIX>   _q;
  FlexVec<double, SM_DFLT_FIX>   _work;

 public:
 DSpline():_nx(0),_k(0),_xb(0.0),_xe(0.0),_inbv(1){}
  ~DSpline(){}

  int   GetNKnots() const { return _k; }

  double GetXb() const { return _xb; }
  double GetXe() const { return _xe; }

  void Init(int npt, int k);
  int Construct(double *xx, double *yy);
  double EvalFun(double x);
  double EvalDer(double x);
  int EvalFunAndDer(double x, double &y, double &dy);

  // might want to have functions which eval an array */

};

#endif

// Local Variables:
// mode: c++
// End:

