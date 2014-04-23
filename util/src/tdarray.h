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

/*
  3D array to store results
 */

#ifndef _TDARRAY_H
#define _TDARRAY_H

#include <new>
#include <assert.h>
#include <string.h>

#include "base.h"

/*
  3D array, width and height are fixed, but the length is automatically adjusted
  
  Each width x height storage unit is a "frame"

  used by the waveform storage

  storage: dimension 1 = height
           dimension 2 = width
           dimension 3 = length
*/

#define MAXLEN 128

template<class T>
class ThreeDArray {
protected:
  int      _width;
  int      _height;
  int      _len;    // allocated length
  int      _loc;     // actually used length

  T*       _vec;

  int      _grow(int nsz) {
    T* nvec;
    try {
      nvec = new T[_width * _height * nsz];
    } catch (std::bad_alloc& ba) {
      return (-1);
    }
    memcpy(nvec, _vec, sizeof(T)*_width*_height*_len);
    delete [] _vec;
    _vec = nvec;
    _len = nsz;
    return 0;
  }

public:
  ThreeDArray():_width(0),_height(0),_len(0),_loc(-1) { _vec=0; } // doomed to cause failure
  ThreeDArray(int w, int h, int l=MAXLEN) {
    assert(w>0 && h>0 && l>0);
    Init(w,h,l);
  }

  // clone the whole thing
  // sizes + content
  ThreeDArray( const ThreeDArray& other) {
    _width = other._width;
    _height = other._height;
    _len    = other._len;
    _loc    = other._loc;
    try {
      _vec = new T[_width * _height * _len ];
    } catch (std::bad_alloc& ba) {
      throw; // hmmm... it would be bad if this happens
    }
    memcpy(_vec, other._vec, sizeof(T)*_width*_height*_len);
  }

  ~ThreeDArray() {
    if (_vec) delete[] _vec;
  }

  /* methods*/
  void Init(int w, int h, int l) {
    if ( w==_width && h==_height && l<_len ) {
      // if same frame size, request less length
      // do nothing, we simply use the previous values
    } else {
      if ( _vec ) delete [] _vec;
      _vec = new T[w*h*l];
      _width = w;
      _height = h;
      _len = l;
    }
    _loc = -1;
  }

  int CurrentLoc() const { return _loc; }

  /* rewind the counter so that we can store new content */
  void Reset() { _loc = -1; }

  /* various methods to store content to the current frame */

  // store one column (height)
  void Store(int h, int w_idx, T* data) {
    assert(h == _height && w_idx < _width);
    for (int jj=0; jj<h; jj++) {
      _vec[_loc*_width*_height + w_idx*_height + jj] = data[jj];
    }
  }
  
  // store one column, but with flexibile starting and ldr to the data vector
  void Store(int h, int w_idx, const int st, const int ldr, T* data) {
    assert(h == _height && w_idx < _width);
    int idx=st;
    for (int jj=0; jj<_height; jj++, idx += ldr) {
      _vec[_loc*_width*_height + w_idx*_height + jj] = data[idx];
    }
  }

  // store the whole frame
  void Store(int frame_len, T* data) {
    assert ( frame_len == _height*_width);
    for (int jj=0; jj<frame_len; jj++) {
      _vec[_loc*_width*_height + jj] = data[jj];
    }
  }

  /* advance the frame counter, and take care of memory allocation */
  int AdvanceFrame() {
    if ( _loc < _len-1 ) return (++_loc);

    // we need to grow allocated storage
    if ( ( _grow( _len*2) ) < 0 ) return (-1);
    return (++_loc);
  }

  /* access methods */

  operator T*() { return _vec; }
  /* not very useful since the caller has to know how data are stored */
  T& operator[](int idx) { return _vec[idx];}

  /* these methods are more useful */

  // return one column at the given frame, and col location
  // returns the size returned
  int GetCol(int len_idx, int w_idx, T* target) const {
    assert( target != NULL );

    if ( len_idx > _loc) return (0); // oops!
    if ( w_idx >= _width) return (0);

    for (int jj=0; jj<_height; jj++) {
      target[jj] = _vec[len_idx*_width*_height + w_idx*_height + jj];
    }
    return (_height);
  }

  // return one row at the given frame, and row location
  int GetRow(int len_idx, int h_idx, T* target) const {
    assert ( target != NULL );

    if ( len_idx > _loc ) return (0);
    if ( h_idx >= _height) return (0);

    for (int jj=0; jj<_width; jj++) {
      target[jj] = _vec[ len_idx*_width*_height + jj*_height + h_idx ];

    }
    return (_width);
  }
  
  // return the (row,col) value across all frames
  int GetSlice(int w_idx, int h_idx, T* target) const {
    assert ( target != NULL );
    
    if ( w_idx >= _width ) return (0);
    if ( h_idx >= _height ) return (0);

    for (int jj=0; jj<=_loc; jj++) {
      target[jj] = _vec[jj*_width*_height + w_idx*_height + h_idx];
    }
    return (_loc+1);
  }
};
#endif

// Local variables:
// mode: c++
// End:
