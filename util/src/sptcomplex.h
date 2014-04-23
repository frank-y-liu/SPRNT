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
  History
  09/01/11    FYL   Split from base.h
*************************************************************************/

#ifndef	_SPT_COMPLEX_H
#define	_SPT_COMPLEX_H

struct  dComplex {
  double re;
  double im;
  double real() { return re; }
  double imag() { return im; }
  dComplex() {};
  dComplex(const double r) : re(r), im(0.0) {};
  dComplex(const double r, const double i) : re(r), im(i) {};
  int operator==(const dComplex& a) { if (a.re == re) return (a.im == im); else return 0; }
  dComplex operator*= (const dComplex& a) {
    double tmp = re;
    re = a.re*re-a.im*im;
    im = a.re*im+a.im*tmp;
    return *this;
  }
  dComplex operator*= (const double& a) {
    re *= a;
    im *= a;
    return *this;
  }
  dComplex operator+= (const dComplex& a) {
    re += a.re;
    im += a.im;
    return *this;
  }
  dComplex operator-= (const dComplex& a) {
    re -= a.re;
    im -= a.im;
    return *this;
  }
};

inline double real( dComplex c ) { return c.real(); }
inline double imag( dComplex c ) { return c.imag(); }

inline dComplex operator* ( dComplex& a,  dComplex& b ) { 
  dComplex c; 
  c.re = a.re*b.re-a.im*b.im;
  c.im = a.re*b.im+a.im*b.re;
  return c;
}

inline dComplex operator* ( dComplex& a,  double b ) { 
  dComplex c; 
  c.re = a.re*b;
  c.im = a.im*b;
  return c;
}

inline dComplex operator* ( double a,  dComplex& b ) { 
  dComplex c; 
  c.re = a*b.re;
  c.im = a*b.im;
  return c;
}

inline dComplex operator/ ( dComplex& a,  dComplex& b ) { 
  dComplex c; 
  double t = b.re*b.re + b.im*b.im;
  c.re = (a.re*b.re+a.im*b.im)/t;
  c.im = (-a.re*b.im+a.im*b.re)/t;
  return c;
}

inline dComplex operator/ ( dComplex& a,  double b ) { 
  dComplex c; 
  c.re = a.re/b;
  c.im = a.im/b;
  return c;
}

inline dComplex operator/ ( double a,  dComplex& b ) { 
  dComplex c; 
  double t = b.re*b.re+b.im*b.im;
  c.re = a*b.re/t;
  c.im = -a*b.im/t;
  return c;
}

#endif

// Local Variables:
// mode: c++
// End:
