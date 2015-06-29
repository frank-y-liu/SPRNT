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
   build the table look up model
   new version: take the horizonal cut line
************************************************************************/

#ifndef _XSMODEL_H
#define _XSMODEL_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "flexvec.h"
#include "dspline.h"
#include "xsection.h"
#include "xstrans.h"

/************************************************************************
   XY cross section
************************************************************************/
class XY_XSection : public XSection {
 protected:
  // the idea is use troughs to build the following  
  int                        _num_segs;    

  DSpline                    _y_sp;    // splines
  DSpline                    _cen_sp;
  DSpline                    _per_sp;
  DSpline                    _w_sp;

  double                     _top_y;
  double                     _bot_y;
  
 public:
 XY_XSection(int id):XSection(id) {
    _top_y = 0.0;
    _bot_y = 0.0;
  }

  ~XY_XSection() { }
  
  // methods for the real usage
  int Build(int num_pairs, double *xx, double *yy, xs_trans *m);
  int Build(int num, double *aa, double *pp, double *yy, double *ww);
	    
  double GetDepth(double a);
  double GetDepthdA(double a);
  double GetWidth(double a);
  double GetCentroid(double a);
  double GetCentroiddA(double a);
  double GetEqFriction(double a);
  double GetEqFrictiondA(double a);

  double CheckCentroiddA(double a);
  double CheckEqFrictiondA(double a);

  double GetAbyDepth(double h);

  void GetCentroid(double a, double &cen, double &cenda);
  void GetEqFriction(double a, double &ef, double &efda);
  void GetHydroRadius(double a, double &r, double &drda);

  int ReachedBankFull(double a) { return ( a >= _aa_bf ); }
  int ReachedMinimalA(double a) { return ( a <  _aa_min); }
};

#endif
// Local Variables:
// mode: c++
// End:

/* End */
