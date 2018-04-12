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

/*****************************************************************
  defines the independent sources 
***************************************************************/

#ifndef _SOURCES_H
#define _SOURCES_H

#include <stdio.h>
#include "base.h"
#include "flexvec.h"
#include "options.h"

/********************************************************
   independent sources, a.k.a. time series
   provide common interface. theoretically we can have
   other types, e.g., sinsq(), sin(), constant() etc
   currently we only have PWL
********************************************************/
class Source {
 protected:
  int                 _parent; 
  SourceType          _type;
  double              _tstart;
  double              _tend;

 public:
  Source(SourceType tp):_parent(-1) {_type=tp;}
  Source(const Source& rhs) { _parent=rhs._parent;_type = rhs._type; }

  virtual ~Source() {}

  // query methods
  double TStart() const { return _tstart; }
  double TEnd()  const  { return _tend;   }

  SourceType Type() { return _type; }
  void SetParent(int j) { _parent = j; }
  int GetParent() const { return (_parent); }
  virtual double Evaluate(double t) { return 0.0;}
  virtual double NextT(double t) { return (FLT_MAX);}

  // returns whether next t, but also the flag when we just cross the breakpoint
  virtual int    GetNextT(double cur, double &nxt) { return 0;}
};

/********************************************************
   PWL: piece-wise-linear
********************************************************/
class PWL : public Source {
 private:
  GrowVec<double,DRSS_PWL_DFT>    _T;    // time
  GrowVec<double,DRSS_PWL_DFT>    _Y;    // value
  GrowVec<double,DRSS_PWL_DFT>    _S;    // slope
  double                          _thresh;
  int                             _size;
  int                             _p1,_p2;

  // make sure the PWL sequence is sane.
  int _check();


  // linear search to find the bracket
  void _lin_search(double t);

  // to make the copy constructor happy! no checking
  double getT(int j) const { return _T.getVal(j); }
  double getY(int j) const { return _Y.getVal(j); }
  double getS(int j) const { return _S.getVal(j); }

 public:
  PWL() : Source(SR_PWL) {_size=0; _p1=-1; _p2=-1; _tstart=-1.0; _tend=-1.0; _thresh=-1.0;}
  ~PWL() {}
  PWL(const PWL& rhs):Source(SR_PWL) { 
    if ( this != &rhs) {
      _p1 = rhs._p1;
      _p2 = rhs._p2;
      _size = rhs._size;
      for (int jj=0; jj<_size; jj++) {
	_Y[jj] = rhs.getY(jj);
	_S[jj] = rhs.getS(jj);
	_T[jj] = rhs.getT(jj);
      }
    }
  }

  double *T() { return &(_T[0]); }
  double *Y() { return &(_Y[0]); }
  double *S() { return &(_S[0]); }
  int Size() const { return _size; }
  double FinalT() const { return (_T.getVal(_size-1)); } // obvious this is only meaningful after
						 // Finalize()
  double FinalY() const { return (_Y.getVal(_size-1)); } // ditto

  // check to see if the time is sequential, we lend it everybody!
  static int CheckTime(int sz, double *t);

  void Scale(double s) {   // scale Y value, T cannot be touched
    assert( s > 0 );
    int jj;
    for (jj=0; jj<_size; jj++) _Y[jj] *= s;
    for (jj=0; jj<_size-1; jj++) _S[jj] *= s;
  }


  /* Method 1:  deprecated
     Assign(0, t1, x1);
     Assign(1, t2, x2);
     rc = Finalize(size)
     if (rc<0) error!  */
  void Assign(int idx, double t, double y);
  int Finalize(int size); // make sure the house is in order

  /* Method2:
     rc = Assign(int size, double*, double*)
     if (rc<0) error!   
     this method can be used to overwrite or concatenate to the existing time series
     example: t1=[0 1 2 3], y1=[....]
              t2=[4 5 6],   y2=[....]
     it's perfectly fine to do:
        Assign(4, t1,y1);
	Assign(3, t2,y2);
     which has the same net effect as calling:
         Assign(7,[0 1 2 3 4 5 6],[.....]);

     this feature is designed for the API mode.
	 
  */
  int Assign(int sz, double *t, double *y);


  /* if t>TEnd(), the time series value will be assumed to be piece-wise-constant */
  double Evaluate(double t); // evaluation given the time point

  double NextT(double t); // return the next time point in the queue
  int GetNextT(double cur, double &nxt); // better way to find the next time point
};

#endif

// Local Variables:
// mode: c++
// End:

// End
