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
   method implementation 
************************************************************************/

#include <math.h>

#include "xsmodel.h"
#include "options.h"

int XY_XSection::Build(int num, double *xx, double *yy, xs_trans *m) {
  int num_sample;
  double *ap, *yp, *lp, *cp, *pp;

  m->Init(_id, num);
  m->Construct(num, &xx[0], &yy[0], OPT.MinB() );
  if ( (num_sample = m->NumSamples()) <= 0 ) {
    fprintf(stderr,"Failed to build table for x-section %d.\n", _id);
    return (-1);
  }

  ap = m->A();
  yp = m->Y();
  lp = m->L();
  cp = m->C();
  pp = m->P();
  
  const int k=3; // was 3rd order spline

  _y_sp.Init(num_sample,k);
  _cen_sp.Init(num_sample,k);
  _per_sp.Init(num_sample,k);
  _w_sp.Init(num_sample,k);

  _y_sp.Construct(&ap[0], &yp[0]);
  _cen_sp.Construct(&ap[0], &cp[0]);
  _per_sp.Construct(&ap[0], &pp[0]);
  _w_sp.Construct(&ap[0], &lp[0]);

  // determine threshold values related to bankful depth and min area
  _aa_bf = ap[num_sample-1]*0.999;
  _aa_min = m->GetMinSampleA()*1e-3;

  // done!

#ifdef DBGSPLINE
  int jj;
  double st=(_aa_bf-_aa_min)/1000;
  printf("spline sampling\n");
  for (jj=0; jj<1000; jj++) {
    printf("=6= %.5e %.5e %.5e %.5e %.5e\n", jj*st, 
	   _w_sp.EvalFun(jj*st), _per_sp.EvalFun(jj*st), _cen_sp.EvalFun(jj*st), _y_sp.EvalFun(jj*st) );
  }
#endif

  return 0;
}

// given A, return depth
double XY_XSection::GetDepth(double a) {
  double h;

  if ( a > _aa_bf ) {
    h = _y_sp.EvalFun( _aa_bf );
    double h2 = _y_sp.EvalFun( 0.98*_aa_bf);
    h = h + ( a - _aa_bf)*(h-h2)/(0.02*_aa_bf);
  } else if ( a < _aa_min) {
    h = _y_sp.EvalFun(a);
  } else {
    h = _y_sp.EvalFun(a);
  }

  return ( h );
}

double XY_XSection::GetDepthdA(double a) {
  double h;

  if ( a > _aa_bf ) {
    h = _y_sp.EvalDer( _aa_bf );
  } else if ( a < _aa_min) {
    h = _y_sp.EvalDer(a);
  } else {
    h = _y_sp.EvalDer(a);
  }

  return ( h );
}

// given A, return centroid
double XY_XSection::GetCentroid(double a) {
  double h;

  if ( a > _aa_bf ) {
    h = _cen_sp.EvalFun( _aa_bf );
    double h2 = _cen_sp.EvalFun( 0.98*_aa_bf );
    h = h + (a-_aa_bf)*(h-h2)/(0.02*_aa_bf);  // assuming linear exp
  } else  if ( a < _aa_min) {
    h = _cen_sp.EvalFun( a );
  } else {
    h = _cen_sp.EvalFun(a);
  }
  return ( h);
}

// given A, return eq friction
double XY_XSection::GetEqFriction(double a) {
  double p, ef;

  if ( a > _aa_bf ) {  // we only print this error message in this function!
    if ( OPT.DebugLevel() >= 2 ) fprintf(stderr,"XY_XSection (%d) have reached bankfull area %.2e (%.2e)\n", _id, a, _aa_bf);
    double p2 = _per_sp.EvalFun( 0.95*_aa_bf);
    p = _per_sp.EvalFun( _aa_bf );
    p = p + (a-_aa_bf)*(p-p2)/(0.05*_aa_bf);
  } else if ( a < _aa_min) {
    if ( OPT.DebugLevel() >= 2 ) fprintf(stderr,"XY_XSection (%d) has small area %.2e (%.2e)\n", _id, a, _aa_min);
    p = _per_sp.EvalFun( a  );
  } else {
    p = _per_sp.EvalFun(a);
  }

  ef = p*cbrt(p/a)/(a*a);

  return ( ef );
}

// given A, return centroid
double XY_XSection::GetCentroiddA(double a) {
  double h;

  if ( a > _aa_bf ) {
    h = _cen_sp.EvalDer( _aa_bf);
  } else if ( a < _aa_min) {
    h = _cen_sp.EvalDer( a );
  } else {
    h = _cen_sp.EvalDer(a);
  }

  return ( h );
}

// given A, return eq friction
double XY_XSection::GetEqFrictiondA(double a) {
  double p, ef, pda,defda;

  if ( a > _aa_bf ) {
    p = _per_sp.EvalFun( _aa_bf);
    pda = _per_sp.EvalDer( _aa_bf);
  } else  if ( a < _aa_min) {
    p = _per_sp.EvalFun( a );
    pda = _per_sp.EvalDer( a );
  } else {
    p = _per_sp.EvalFun(a);
    pda = _per_sp.EvalDer(a);
  }

  ef = p*cbrt(p/a)/(a*a);
  defda = ef*( 4.0*pda/(3.0*p) - 7.0/(3.0*a) );

  return ( defda );
}


void XY_XSection::GetCentroid(double a, double &cen, double &cenda) {
  double fcen, fdcen;

  if ( a > _aa_bf ) {
    fcen = _cen_sp.EvalFun(_aa_bf);
    double c2 = _cen_sp.EvalFun( 0.95*_aa_bf );
    fcen = fcen + (a-_aa_bf)*(fcen-c2)/(0.05*_aa_bf);
    fdcen = _cen_sp.EvalDer(_aa_bf);
  } else if ( a < _aa_min) {
    fcen = _cen_sp.EvalFun( a );
    fdcen = _cen_sp.EvalDer( a );
  } else {
    fcen = _cen_sp.EvalFun(a);
    fdcen = _cen_sp.EvalDer(a);
  }

  // if ( fcen < 0 && ABS(fcen)<1e-6) fcen = -fcen;
  cen = fcen;
  cenda = fdcen;

}

void XY_XSection::GetEqFriction(double a, double &ef, double &efda) {
  double p, pda, fef, fdef;
  
  if ( a > _aa_bf ) {
    double p2 = _per_sp.EvalFun( 0.95*_aa_bf);
    p = _per_sp.EvalFun(_aa_bf);
    p = p + (a-_aa_bf)*(p-p2)/(0.05*_aa_bf);
    pda = _per_sp.EvalDer(_aa_bf);

  } else if ( a < _aa_min) {
    p = _per_sp.EvalFun( a );
    pda = _per_sp.EvalDer( a );
  } else {
    p = _per_sp.EvalFun(a);
    pda = _per_sp.EvalDer(a);
  }

  fef = p*cbrt(p/a)/(a*a);
  fdef = (4.0/3.0)*fef*pda/p - (7.0/3.0)*fef/a;

  ef = fef;
  efda = fdef;

}


double XY_XSection::CheckCentroiddA(double a) {
  double a2 = a + _aa_min;
  double c1, c2;

  c1 = this->GetCentroid(a);
  c2 = this->GetCentroid(a2);
  return ( (c2-c1)/_aa_min );
  
}

double XY_XSection::GetWidth(double a) {
  double w;

  if ( a > _aa_bf ) {
    w = _w_sp.EvalFun( _aa_bf );
    double w2 = _w_sp.EvalFun( 0.98*_aa_bf );
    w = w + (a-_aa_bf)*(w-w2)/(0.02*_aa_bf);  // assuming linear exp
  } else  if ( a < _aa_min) {
    w = _w_sp.EvalFun( a );
  } else {
    w = _w_sp.EvalFun(a);
  }
  return (w);
}

double XY_XSection::CheckEqFrictiondA(double a) {
  double a2 = a + 1e-3;
  double e1, e2;

  e1 = this->GetEqFriction(a);
  e2 = this->GetEqFriction(a2);

  return ( (e2-e1)/1e-3 );

}

/* given A, calculate hydraulic radius and derivative */
void XY_XSection::GetHydroRadius(double a, double &r, double &drda) {
  double p, dp;

  p = _per_sp.EvalFun( a );
  dp = _per_sp.EvalDer( a );

  r = a/p;
  drda = ( 1 - a *dp/p)/p;
}

// This method uses the coord method, could be expensive
double XY_XSection::GetAbyDepth(double h) {
  double a0, a1, h0, h1, aa, htest;
  int jj;
  const double res = 1e-6;
  const int maxiter = 50; 

  h0 = _y_sp.EvalFun(_aa_min);
  h1 = _y_sp.EvalFun(_aa_bf);
  
  if ( h < h0 ) {         // too small, return min
    aa = _aa_min; 
  } else if ( h > h1 ) {  // too high, return projection
    double h2 = _y_sp.EvalFun(_aa_bf*0.98);
    aa = _aa_bf + ( 0.02*_aa_bf *(h-h2)/(h1-h2));
  } else {                // otherwise do a binary search
    a0 = _aa_min;
    a1 = _aa_bf;
    for (jj=0; jj<maxiter; jj++) {
      // will return crappy results if the iteration doesn't converge within the maxiter 
      aa = a0 + (a1-a0)*(h-h0)/(h1-h0);
      htest = _y_sp.EvalFun(aa);
      
      if ( fabs(htest-h)< res ) break;
      if ( htest < h ) {
	a0 = aa;
	h0 = htest;
      } else if ( htest > h ) {
	a1 = aa;
	h1 = htest;
      }
    }
  }
  return (aa);
}

// End


