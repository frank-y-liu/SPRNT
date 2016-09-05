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

/* header files for the simple parser */

#ifndef _SIM_PARSE_H
#define _SIM_PARSE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

// max number of characters per line
#define MAX_LINE_LENGTH 512
#define MAX_WORD_LENGTH 64
#define SHORT_ARRY_LENGTH 256

// The sequence here has to be exacly same as Spt_Descriptor!
// except the added values
typedef enum {
  p_options = 0,
  p_node,
  p_segment,
  p_qsource,
  p_boundarycondition,
  p_lateralsource,
  p_junction,
  def_not_found,
  def_error,
  def_block            /* add this later */
} FirstDef;

typedef enum {
  trapezoidal = 1,
  rectangular,
  xy,
  timeseries
} SecondDef;

typedef enum {
  ft_ascii,
  ft_integer,
  ft_real,
  ft_block
} FType;

typedef enum {
  no = 0,
  yes,
  secondary
} Required;

typedef struct {
  const char      *_key;
  const Required  _required;
  const FType     _type;
} Descriptor;

#if 0
// simple memory manager for double array to avoid repeated free/malloc
// The user is reponsible of keeking track of the indices
class SimpleDblArray {
 private: 
  double  *_xx;
  double  *_yy;
  double   _scale;
  int      _alloc_size;
  int      _size;
  
  inline void _grow(int sz) {   // check the size and grow if needed
    if ( sz < _alloc_size ) return;
    double *tmpx, *tmpy;
    int newsz = 2*sz;
    tmpx = (double*)malloc(sizeof(double)*newsz);
    tmpy = (double*)malloc(sizeof(double)*newsz);
    memcpy(tmpx, _xx, sizeof(double)*_alloc_size);
    memcpy(tmpy, _yy, sizeof(double)*_alloc_size);
    free(_xx); free(_yy);
    _xx = tmpx;
    _yy = tmpy;
    _alloc_size = newsz;
  }

 public:
  SimpleDblArray(int sz=512) {
    _alloc_size = sz;
    _xx = (double*) malloc(sizeof(double)*sz);
    _yy = (double*) malloc(sizeof(double)*sz);
    _size = 0;
    _scale = 1.0;
  };

  ~SimpleDblArray() {
    if (_xx) free(_xx);
    if (_yy) free(_yy);
  }

  double *X() { return _xx; }
  double *Y() { return _yy; }
  double& ithX(int i) { _grow(i); return _xx[i]; }
  double& ithY(int i) { _grow(i); return _yy[i]; }
  void Reset() { _size = 0; _scale = 1.0; }
  const int Size() const { return _size; }
  void SetScale(double s) { _scale = s; }
  const double Scale() const { return _scale; }

};
#else
/*
  Simple N-tuple, by default N=2
  default size = 512
  This template is intended to be used to temporary store data, therefore it is not bullet
  proof. Avoid letting non-friends to use it

  To assigned the data: use ithVal() method
  getIthArray() is provided to quickly access the data, but it is not safe
  
  Reset() will reset the template for next use. Memory is not freed till the instance is
  deleted

 */
template<class T, unsigned int NSEQ=2>
class simpleTuple {
 protected:
 T     **_data;
 T       _scale;

 int     _alloc_size;
 int     _size;

 inline void _grow(int sz) {
   if ( sz < _alloc_size) return;
   T *tmp;
   int newsz = 2*sz; // growth policy
   for (unsigned int jj=0; jj<NSEQ; jj++) {
     tmp = (T*)malloc(sizeof(T)*newsz);
     memcpy(tmp, _data[jj], sizeof(T)*_alloc_size);
     free(_data[jj]);
     _data[jj] = tmp;
   }
   _alloc_size = newsz;
 }

 public:
 simpleTuple(int sz=512) {
   _data = (T**)malloc(NSEQ*sizeof(T*));
   for (unsigned int jj=0; jj<NSEQ; jj++)  _data[jj] = (T*)malloc(sizeof(T)*sz);
   _alloc_size = sz;
   _size = 0;
   _scale = (T)1.0;
 }

 ~simpleTuple() {
   for (unsigned int jj=0; jj<NSEQ; jj++)  if (_data[jj]) free(_data[jj]);
   free (_data);
 }

 // methods:
  T &ithVal(int jj, int kk) { assert((unsigned int)jj<NSEQ); _grow(kk); return (_data[jj][kk]); }
  T *getIthArray(int jj) { assert((unsigned int)jj<NSEQ); return (_data[jj]); } /* be careful */

 void Reset() { _size = 0; _scale=(T)1.0; }
 const int Size() const { return _size; }
 void SetScale(T s) { _scale = s; }
 const T Scale() const { return _scale; }

};

#endif

// class to store the node names, by the given index
class NameStore {
 private:
  char   **_store;
  int      _allocated_size;
  int      _word_length;
  int      _used;

 public:
  NameStore() {
    _store = NULL;
    _allocated_size = -1;
    _word_length = -1;
    _used = -1;
  }

  ~NameStore() {
    if (_store) {
      for (int jj=0; jj<_allocated_size; jj++) if (_store[jj]) free(_store[jj]);
      free(_store);
    }
  }

  // simple methods, Init(), Insert() and Get()
  // if you call the function before Init(), you will get a seg fault
  void Init(int sz, int word_length=MAX_WORD_LENGTH) {
    assert(sz>0 && word_length>0);
    _used = -1;
    _allocated_size = sz+10; // add some extra space
    _word_length = word_length;
    _store = (char**)malloc(sizeof(char*)*_allocated_size);
    for (int jj=0; jj<_allocated_size; jj++) _store[jj]=(char*)malloc(sizeof(char)*word_length);
  }

  void Insert(char *c, int w) {
    assert( w>=0 && w<_allocated_size);
    strncpy(_store[w], c, _word_length);
    _used = _used > w ? _used : w;  // take the max
  }
  char* Get(int w) {
    assert( w>=0 && w<_allocated_size);
    return _store[w];
  }

  const int MaxIndex() const { return _used; }
  char** Store() { return _store; }

};

// Sring buffer to store the results of parsing
class StrBuffer {
 private:
  char    *_buffer;
  char    *_current;
  char     _header[MAX_LINE_LENGTH];
  char    **_keys;
  char    **_values;
  int      _pair_size;
  int      _pair_used;
  int      _size;
  int      _used;
  int      _line_number;

  void _check(int sz);  // check the memory, and grow if needed

  void _free_pairs();   // free allocated memory

  void _allocate_pairs(int sz); // allocate memory

  void _check_pair(int sz) {
    if ( sz < _pair_size) return;
    sz = 2*sz;
    _free_pairs();
    _allocate_pairs(sz);
  }

 public:
  StrBuffer(int init_size, unsigned int num_pairs=2) { 
    assert(init_size>0);
    _buffer = (char*) malloc(sizeof(char)* init_size);
    _size = init_size;
    _allocate_pairs(num_pairs);
    Reset();
  }

  ~StrBuffer() {
    if ( _buffer ) free(_buffer);
    _buffer = 0; 
    _size = 0;
    _free_pairs();
  }

  /* Typical usage should be:
     Reset(), Ready(), Copy_add_Space(), Finish(), and repeat 
     use HasContent() to see if there is anything
  */
  inline void Reset() { _used=0; _line_number=0; _current = NULL; _buffer[0]=0;}
  inline void Ready(int ln) { _current = _buffer; _line_number = ln; }
  void Copy_add_Space(char *s);
  inline void Finish() { if (_current) {*_current++ = 0; _used++;} }
  
  /* utility functions, HasContent() tells if there is anything inside */  
  const int LineNumber() const { return _line_number; }
  const int NumPairs() const { return _pair_used; }
  inline int HasContent() { return (_current!=NULL); }
  inline void Dump(FILE *out) { if (out && HasContent()) fprintf(out, "%d: %s\n",_line_number,_buffer); }

  inline int Cmp_Head(const char *t2) { // compare to see if the header matches the given char
    int rc = 1;
    char *t1 = _header;
    while (1) {
      if (*t2==0 || rc==0 ) break;
      rc *= (toupper(*t1++) == *t2++ ? 1 : 0); /* logic AND, need all characters */
    }
    return (rc==1);
  }

  char* Get_ith_Value(const char *k, int j);
  char* Find_Value(const char *k);
  int Separate(FILE *out);  // parsing the results and put them into the pairs in
			    // key-value pair
  
};


// for convenience
typedef simpleTuple<double>   SimpleDblArray;
typedef simpleTuple<double,4> SimpleDblQuad;

#endif

// Local Variables:
// mode: c++
// End:
