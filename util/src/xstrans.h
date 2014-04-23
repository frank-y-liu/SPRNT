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
   method to translate arbitrary cross section into sampled 
   tables. Apply spline fitting
************************************************************************/

#ifndef _XSTRANS_H
#define _XSTRANS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "flexvec.h"

// for memory allocation only
const int XSDFLT_FIX = 255;

// controls the behavior of the code
const double                ANG_THRESH = 0.001;
const double                UNIT = 0.001; // default 0.001
const int                   YGRID = 10; // default 10
const int                   SGRID = 10; // default 10

typedef struct _segment {
  long int _x1;
  long int _y1;
  long int _x2;
  long int _y2;
} segment;

/* holds the method to convert XY x-section to spline-fitted */
class xs_trans {
 protected:
  int                           _name;
  FlexVec<double, XSDFLT_FIX>   _my_aa;
  FlexVec<double, XSDFLT_FIX>   _my_ll;
  FlexVec<double, XSDFLT_FIX>   _my_pp;
  FlexVec<double, XSDFLT_FIX>   _my_yy;
  FlexVec<double, XSDFLT_FIX>   _my_cc;

  FlexVec<double, XSDFLT_FIX>   _aa;
  FlexVec<double, XSDFLT_FIX>   _ll;
  FlexVec<double, XSDFLT_FIX>   _pp;
  FlexVec<double, XSDFLT_FIX>   _cc;
  FlexVec<double, XSDFLT_FIX>   _yy;

  FlexVec<long int, XSDFLT_FIX>  _xxin;
  FlexVec<long int, XSDFLT_FIX>  _yyin;
  FlexVec<long int, XSDFLT_FIX>  _rxx;
  FlexVec<long int, XSDFLT_FIX>  _ryy;
  FlexVec<int, XSDFLT_FIX>       _tmp;
  
  FlexVec<long int, XSDFLT_FIX> _working;
  FlexVec<long int, XSDFLT_FIX> _prev;

  FlexVec<segment, XSDFLT_FIX>  _segs;

  long int                    _min_y;
  long int                    _max_y;

  int                         _num_segs;
  int                         _num_pts;
  int                         _num_samples;
  int                         _final_samples;

  // convert to long int
  int _record_into_int(int num, double *xin, double *yin, double scale); 
  // cone method to reduce size
  int _cmethod(double threshold);
  // facility: calculate angle
  double _angle(long int x1, long int y1, long int xs, long int ys); 
  // insert into segments
  int _insert2segs(); 
  // find intersection with given yy, results in working, return number
  int _find_intersects(long int yy, long int *working);
  // get len, also does scaling
  double _get_len(int num_int, long int *working, double unit); 
  // get perimeter, assuming trapezoidal, also does scaling
  double _get_pp(int num_now, long int *now, int num_prev, long int *prev, double dy, double unit);
  // find min and max, stored them
  void _find_min_max();
  // initialize internal data storage
  void _init_internal();
  // a quick check to make sure the integrity of the channel is good
  int _quick_check();

 public:
   xs_trans():_num_segs(0),_num_pts(0),_num_samples(0),_final_samples(0) {}
  ~xs_trans() {}

  // initialization
  void Init(int cname, int sz);
  // do real work, once it is done, results are in arrays
  int Construct(int num, double *x, double *y, double bot_thresh);
  
  // access method, use them to get details of the x-section after calling Construct
  int NumSamples() const { return _final_samples; }
  double GetMinSampleA() { return _my_aa[1]; }
  double *A() { return &_aa[0]; }
  double *Y() { return &_yy[0]; }
  double *L() { return &_ll[0]; }
  double *C() { return &_cc[0]; }
  double *P() { return &_pp[0]; }

  // These methods are for debugging only
  void PrintDenseAll(FILE *F, const char *h="-D-");
  void PrintSparseAll(FILE *F, const char *h="-S-");

};

#endif

// Local Variables:
// mode: c++
// End:

/* end */
