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

/****************************************************************************
  cross section definition
  09/01/11    so far for R_ and T_ types only
****************************************************************************/

#include <math.h>
#include <stdio.h>
#include "xsection.h"
#include "options.h"

// R type, rectangular
inline double R_XSection::GetDepth(double a) {
  return a/_b0;
}

inline double R_XSection::GetDepthdA(double a) {
  return 1.0/_b0;
}

inline double R_XSection::GetWidth(double a) {
  return _b0;
}

inline double R_XSection::GetCentroid(double a) {
  return 0.5*a*a/_b0;
}

inline double R_XSection::GetCentroiddA(double a) {
  return a/_b0;
}

inline double R_XSection::GetEqFriction(double a) {
  double t = _b0 + 2.0*a/_b0;
  return (  t*cbrt(t/a)/(a*a) );
}

inline double R_XSection::GetEqFrictiondA(double a) {
  double pp = _b0+2.0*a/_b0;
  double ef = pp*cbrt(pp/a)/( a*a );
  double def = ef*( 8.0/(3.0*_b0*pp) - 7.0/(3.0*a) );
  return def;
}

void R_XSection::GetCentroid(double a, double &cen, double& cenda) {
  cenda = a/_b0;
  cen = 0.5*a*cenda;
}

void R_XSection::GetEqFriction(double a, double &ef, double &efda) {
  double pp = _b0+2.0*a/_b0;
  ef = pp*cbrt(pp/a)/(a*a);
  efda = ef*( 8.0/(3.0*_b0*pp) - 7.0/(3.0*a) );
}

void R_XSection::GetHydroRadius(double a, double &r, double &drda) {
  double p;   // wetted perimeter
  p = 2*a/_b0 + _b0;
  r = a/p;
  drda = (1.0 - a/p *( 2/_b0))/p;
}

// T type, trapezoidal 
inline double T_XSection::GetDepth(double a) {
  return ( sqrt(_effb0*_effb0 + a/_s) - _effb0 );
}

inline double T_XSection::GetDepthdA(double a) {
  return ( 1/( 2 * _s * sqrt( _effb0*_effb0 + a/_s) ) );
}

inline double T_XSection::GetWidth(double a) {
  return ( _b0 + 2*sqrt(_effb0*_effb0 + a/_s) - 2*_effb0 );
}

inline double T_XSection::GetCentroid(double a) {
  double dd = sqrt(_effb0*_effb0 + a/_s) - _effb0;
  return dd*dd*( _b0/2.0 + dd*_s/3.0);
}

inline double T_XSection::GetCentroiddA(double a) {
  double sqdd = sqrt(_effb0*_effb0 + a/_s);
  double dcda = ( (sqdd - _effb0)*( _effb0 + 0.5*(sqdd-_effb0) ) /sqdd );
  return dcda;
}

inline double T_XSection::GetEqFriction(double a) {
  double dd = sqrt( _effb0*_effb0 + a/_s) - _effb0;
  double tt = 2*sqrt(1+_s*_s)*dd + _b0;
  double ef = (tt*cbrt(tt/a))/(a*a);
  return ef;
}

double T_XSection::GetEqFrictiondA(double a) {
  double sqdd = sqrt( _effb0*_effb0 + a/_s);
  double dd = sqdd - _effb0;
  double sp1q = sqrt(1+_s*_s);
  double tt = (2*sp1q*dd + _b0);
  double ef = tt*cbrt(tt/a)/(a*a);
  double efda = ef*( (4.0/3.0)*sp1q/(_s*sqdd*(2*sp1q*dd+_b0)) - (7.0/3.0)/a );

  return efda;
}

void T_XSection::GetCentroid(double a, double &cen, double& cenda) {
  double sqdd = sqrt( _effb0*_effb0 + a/_s );
  double dd = sqdd - _effb0;
  cen = dd*dd*(_b0/2.0 + dd*_s/3.0);
  cenda = dd*( _effb0 + 0.5*dd ) / sqdd;
}

void T_XSection::GetEqFriction(double a, double &ef, double &efda) {
  double sqdd = sqrt( _effb0*_effb0 + a/_s);
  double dd = sqdd - _effb0;
  double sp1q = sqrt(1+_s*_s);
  double tt = (2*sp1q*dd + _b0);
  ef = tt*cbrt(tt/a)/(a*a);
  efda = ef*( (4.0/3.0)*sp1q/(_s*sqdd*(2*sp1q*dd+_b0)) - (7.0/3.0)/a );
}

void T_XSection::GetHydroRadius(double a, double &r, double &drda) {
  double t1 = sqrt(_effb0*_effb0 + a/_s);
  double t2 = sqrt(1+_s*+_s);
  r = _b0 + 2 * t2 * ( t1 - _effb0);
  drda = t2/(_s*t1);
}

// Local Variables:
// mode: c++
// End:


/* end */
