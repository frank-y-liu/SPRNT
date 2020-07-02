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

#ifndef _FLEXVEC_H
#define _FLEXVEC_H

#include <assert.h>
#include <memory.h>
#include "base.h"

/* 
   Vector with flexible length
   resize() will change the size, and wipe out the existing content
   grow() will grow, but also keep the existing content
   Caveat: frequent small increase may cause performance degradation, use GrowVec instead
*/

template<class T, int FIXSZ=MAXFIXSZ> 
class FlexVec {
 protected:
  T   _fixmem[FIXSZ];
  int _asz;
  T*  _vec;
 //if not large enough allocates more than requested and recopies the old data into new memory
  int do_grow( int sz ) { 
    T* newvec = new T[sz + _asz]; // at least doubles the allocated size
    for( int k=0; k!=_asz; ++k) newvec[k] = _vec[k];
    if (_vec != (T*)_fixmem) delete[] _vec;
    _vec = newvec;
    _asz += sz;                            
    return _asz;
 }

 public:
  /* constructors and destructor */
 FlexVec() : _asz(FIXSZ), _vec(_fixmem) {};
 FlexVec( int sz ) { 
    if (sz <= FIXSZ) {
      _vec = _fixmem;
      _asz = FIXSZ;
    } else {
      _vec = new T[sz];
      _asz = sz;
    }
 };
 FlexVec( const FlexVec<T,FIXSZ>& ) { assert(0); } // disallow copy constructor
 ~FlexVec() {
   if (_vec != (T*)_fixmem) delete[] _vec; 
 };

 /* access methods */
 operator T*() { return _vec; }
 T& operator[](int i) { assert(i<_asz); assert(i>=0); return _vec[i]; }; 

 int size( int sz ) { // if not large enough allocates exactly as requested
   if ( sz <= _asz ) return _asz; 
   if (_vec != (T*)_fixmem) delete[] _vec;
   _vec = new T[sz];                
   _asz = sz;                            
   return _asz;
 }

 int resizep( int sz ) { // if not large enough allocates more than requested
   if ( sz <= _asz ) return _asz; 
   if (_vec != (T*)_fixmem) delete[] _vec;
   _vec = new T[sz + _asz]; // at least doubles the allocated size                
   _asz += sz;                            
   return _asz;
 }

 inline int grow (int sz) { return (sz<_asz)?_asz:do_grow(sz); }

};

/*
  Flexible vector with automatic growth strategy, derived from FlexVec
 */
template<class T, int FIXSZ=MAXFIXSZ>
class GrowVec : public FlexVec<T,FIXSZ> {
  public:
  GrowVec()         : FlexVec<T,FIXSZ>()   {};
  GrowVec( int sz ) : FlexVec<T>(sz) {};
  GrowVec( const GrowVec<T,FIXSZ>& ) { }
  T& operator[](int i) { 
    FlexVec<T,FIXSZ>::grow(i);
    return FlexVec<T,FIXSZ>::_vec[i];
  }
  T getVal(int i) const { 
    return FlexVec<T,FIXSZ>::_vec[i];
  }
};

/*
  GrowVecPair, a pair of GrowVec
 */
template<class T, int FIXSZ=MAXFIXSZ>
class GrowVecPair {
protected:
  GrowVec<T, FIXSZ>  _first;
  GrowVec<T, FIXSZ>  _second;
  
public:
  GrowVecPair()      :_first(), _second() {}
  GrowVecPair(int sz):_first(sz),_second(sz) {}
  GrowVecPair( const GrowVecPair<T,FIXSZ>& ) {}
  T& First(int i) { return ( _first[i] );  }
  T& Second(int i) { return ( _second[i] ); }
  T getFirstVal(int i) const { return _first.getVal(i); }
  T getSecondVal(int i) const { return _second.getVal(i); }
};

typedef GrowVecPair<double>  doublePair;

#endif

// Local Variables:
// mode: c++
// End:

