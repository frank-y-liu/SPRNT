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

/**************************************************************************
  implementations for source
**************************************************************************/
#include <assert.h>
#include "sources.h"

int PWL::CheckTime(int sz, double *t) {
  double thresh = (t[sz-1] - t[0])*OPT.Epsilon();
  if ( (t[sz-1]-t[0])<OPT.Epsilon() ) thresh = OPT.Epsilon();
  if ( sz < 2 ) return (-1);
  
  for (int jj=1; jj<sz; jj++) {
    if ( t[jj] - t[jj-1] < thresh ) return (-1);
  }
  return 0;
}

// make sure the time sequence is monotonically increasing, also calculate the slopes
int PWL::_check() { 

  _thresh = (_T[_size-1] - _T[0])*OPT.Epsilon();
  
  // initialize to the first bracket
  _p1 = 0;
  _p2 = 1;
  
  if ( _T[0] < 0.0 || _thresh < 0.0 ) { // doodoo, we don't like it
    _size = 0;
    return -1;
  }
  
  if ( PWL::CheckTime(_size, &_T[0]) < 0 ) {
    _size = 0;
    return (-1);
  }

  // we are good to go
  for (int i=1; i<_size; i++) {       // walk through to compute slope 
    _S[i-1] = (_Y[i]-_Y[i-1])/(_T[i]-_T[i-1]);
  }

  _tstart = _T[0];
  _tend = _T[_size-1];
  
  return 0;
}

/* the purpose of _p1 and _p2 is to speed up the evaluation speed
   since the evaluation time is most likely sequential, we store the two previous "knots"
   so that we don't have to search for the "knots" every time.
*/ 
void PWL::_lin_search(double t) {
  if ( t < _T[0] ) {
    _p1 = 0;
    _p2 = 1;
    return;
  }
  if ( t>_T[_size-1] ) {
    _p1 = _size-2;
    _p2 = _size-1;
    return;
  }
  while ( t > _T[_p2] ) {
    _p2++;
    _p1++;
  }
}

// as the name says
double PWL::Evaluate(double t) {
  double p=0.0;
  
  if ( !(_T[_p1]<t && _T[_p2]>t) ) _lin_search(t);  // and modify _p1, _p2
  
  if ( t>_T[_size-1] ) {
    p = _Y[_size - 1];     // use the last value if t is too large
  } else  {
    p = _Y[_p1] + _S[_p1]*(t-_T[_p1]);
  }
  if (OPT.DebugLevel() >= 3 ) {
    printf("    PWL %lx, Y(%.2e)=%.2e\n", (unsigned long)this, t, p );
  }
  return p;
}

double PWL::NextT(double t) {
  if ( !(_T[_p1]<t && _T[_p2]>t) ) _lin_search(t);
  if ( _T[_p2] - t < TUNIT ) {  // if we are within striking limit of the next point, we
    // move over
    if ( _p2 == _size-1 ) return (FLT_MAX); 
    return (_T[_p2+1]);
  }
  return (_T[_p2]);
}

// better way to find the next time point
int PWL::GetNextT(double cur, double &nxt) {
  int rc;

  if ( !(_T[_p1]<cur && _T[_p2]>cur) ) _lin_search(cur);
  
  nxt = FLT_MAX;
  rc = 0;
  if ( _T[_p2] - cur < TUNIT ) {  
    // if we are within striking limit of the next point, move over
    if ( _p2 == _size-1 ) { nxt = FLT_MAX;} // we hit the last one
    else { rc = 1; nxt=_T[_p2+1]; }
  } else {
    rc = 0;
    nxt = _T[_p2];
  }
  return rc;
}

// return 0 means everything is fine
// if there is time overlap, the new value will overwrite the existing time series
// if there is no overlap, the new series is concatenated to the existing series
int PWL::Assign(int sz, double *t, double *y) {
  int jj;
  if ( sz <= 0 ) { return -1;}
  if ( _size == 0 ) {
    for (jj=0; jj<sz; jj++) {
      _T[jj] = t[jj];
      _Y[jj] = y[jj];
    }
    _size = sz;
  } else {
    int cntr;
    for (cntr=0; cntr <_size; cntr++) {
      if ( (_T[cntr]-t[0]) > _thresh ) break;
    }
    _size -= (_size - cntr);
    for (jj=0; jj<sz; jj++) {
      _T[jj+cntr] = t[jj];
      _Y[jj+cntr] = y[jj];
      _size++;
    }
  }

  return ( _check() );
}

// this method is less flexible than the other overloaded version
void PWL::Assign(int idx, double t, double y) {
  assert(idx >= 0);
  _T[idx] = t;
  _Y[idx] = y;
}

int PWL::Finalize(int size) {
  assert(size>0);
  _size = size;
  return ( _check() );
}

//End
