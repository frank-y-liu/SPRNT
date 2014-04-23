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

/*********************************************************************
 interface to BLAS routines

 make sure to define LINUX 
*********************************************************************/

#ifndef _BLASPR_H
#define _BLASPR_H

#include "base.h"
#include "sptcomplex.h"

#ifdef LINUX
#define FORTRAN(a)  a##_
#else
#define FORTRAN(a)  a
#endif

extern "C" {

double FORTRAN(d1mach)(const int*);
void   FORTRAN(dcopy)(int* n, const double* dx, const int* incx, double* dy, const int* incy);
double FORTRAN(dnrm2)(int*, const double*, const int*);
void   FORTRAN(dscal)(int*, const double*, const double*, const int*);
void   FORTRAN(zscal)(int*, const dComplex*, const dComplex*, const int*);
double FORTRAN(ddot) (int*, const double*, const int*, const double*, const int*);
void   FORTRAN(daxpy)(int*, const double*, const double*, const int*, double*, const int*);
int    FORTRAN(idamax)(int*, const double*, const int*);
void   FORTRAN(zscal)(int*, const dComplex*, const dComplex*, const int*);

}

inline void Sscal( int n, const double alp, double* vec, int inc=1) { 
  FORTRAN(dscal)(&n, &alp, vec, &inc); 
}
/* BLAS routine sscal
*  Purpose
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to 1.
*  ===================================================================== */

inline void Sscal( int n, const dComplex alp, dComplex* vec, int inc=1) { 
  FORTRAN(zscal)(&n, &alp, vec, &inc); 
}
/* BLAS routine sscal
*  Purpose
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to 1.
*  ===================================================================== */


#endif // BLASPR_H

// Local Variables:
// mode: c++
// End:

