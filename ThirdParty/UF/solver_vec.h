/*********************************************************************
  Copyright (C) 2011, 2014 International Business Machines
  All Rights Reserved

  Author:  Frank Liu, IBM
*********************************************************************/

#ifndef _SOLVER_VEC_H
#define _SOLVER_VEC_H

#include <assert.h>

/****** 
   vector with flexible length
   size() and resize() will change the size, and wipe out the existing content
   grow() will grow, but also keep the existing content
   no need to worry about memory 
*******/

#define MAXFIXEDSIZE 512

template<class T, int fixsize=MAXFIXEDSIZE> 
class flexvec {
protected:
  T   _fixmem[fixsize];
  int _asz;
  T*  _vec;

  // _do_grow() will grow and coppy the old data
  int _do_grow( int sz ) { 
    T* newvec = new T[sz + _asz]; // at least doubles the size
    for( int k=0; k!=_asz; ++k) newvec[k] = _vec[k];
    if (_vec != (T*)(&_fixmem[0])) delete[] _vec;
    _vec = newvec;
    _asz += sz;                            
    return _asz;
  }
  
public:
  flexvec():_asz(fixsize), _vec(&_fixmem[0]) {};
  flexvec(int sz) { 
    if (sz <= fixsize) {
      _vec = &_fixmem[0];
      _asz = fixsize;
    } else {
      _vec = new T[sz];
      _asz = sz;
    }
  };
  
  // disallow any copy constructor
  flexvec( const flexvec<T,fixsize>& ) { assert(0); }
  ~flexvec() {
    if (_vec != (T*)(&_fixmem[0]) ) delete[] _vec; 
  };
  
  operator T*() { return _vec; }
  T& operator[](int i) { assert(i<_asz); assert(i>=0); return _vec[i]; }; 

  // size() and resize() should be used before any content is in place
  int size( int sz ) { 
    if ( sz <= _asz ) return _asz; 
    if (_vec != (T*)(&_fixmem[0])) delete[] _vec;
    _vec = new T[sz];                
    _asz = sz;                            
    return _asz;
  }
  
  int resize( int sz ) { 
    if ( sz <= _asz ) return _asz; 
    if (_vec != (T*)(&_fixmem[0])) delete[] _vec;
    _vec = new T[sz + _asz]; 
    _asz += sz;                            
    return _asz;
  }
  
  // grow is safe to use after the content has been added
  inline int grow (int sz) {
    return (sz<_asz)?_asz:_do_grow(sz); 
  }
};

#endif

// Local Variables:
// mode: c++
// End:

