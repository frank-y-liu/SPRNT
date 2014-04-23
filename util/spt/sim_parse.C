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

// methods for parsing netlists

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "sim_parse.h"

// methods for StrBuffer
void StrBuffer::_check(int sz) {
  if (_used + sz < _size - 9 ) return;
  int sz1 = _used + sz + 9;
  int sz2 = 2*_size;
  sz1 = sz1 >= sz2? sz1 : sz2;  // find the max
  char *tmp = (char*) malloc( sizeof(char)*sz1 );
  strncpy(tmp, _buffer, _size);
  free(_buffer);
  _buffer = tmp;
  _current = _buffer;
  _current += _used;   // reassign the pointer and also move it
  _size = sz1;
}

void StrBuffer::_free_pairs() {
  int jj;
  for (jj=0; jj<_pair_size; jj++) {
    if (_keys[jj]) free(_keys[jj]);
    if (_values[jj]) free(_values[jj]);
  }
  if (_keys) free(_keys);
  if (_values) free(_values);
  _keys = NULL;
  _values = NULL;
}

void StrBuffer::_allocate_pairs(int sz) {
  assert(sz>0);
  _keys = (char**)malloc(sizeof(char*)*sz);
  _values = (char**)malloc(sizeof(char*)*sz);
  for (int jj=0; jj<sz; jj++) {
    _keys[jj] = (char*)malloc(sizeof(char)*MAX_LINE_LENGTH);
    _values[jj] = (char*)malloc(sizeof(char)*MAX_LINE_LENGTH);
  }
  _pair_size = sz;
}

void StrBuffer::Copy_add_Space(char *s) {
  int sz=0;
  char *tmp;
  
  sz = strlen(s);
  _check(sz);
  
  tmp = s;
  while (1) {
    if (*tmp==0) break;
    *_current++ = *tmp++;
  }
  *_current++ = ' ';
  _used += sz+1;  /* extra 1 for space */
};

char* StrBuffer::Get_ith_Value(const char *k, int j) {
  const char *t1, *t2;
  int rc=1;
  
  if ( j < 0 || j > _pair_used-1 ) return NULL;
  
  t1 = _keys[j];
  t2 = k;
  while (1) {
    if (*t2==0 || rc==0 ) break;
    rc *= (toupper(*t1++) == *t2++ ? 1 : 0);
  }
  
  if (rc ==1) return (_values[j]);
  else return (NULL);
}

char* StrBuffer::Find_Value(const char *k) {
  int jj,rc;
  char *t1;
  const char *t2;
  
  if (_pair_used==0) return NULL;
  for (jj=0; jj<_pair_used; jj++) {
    rc = 1;
    t1 = _keys[jj];
    t2 = k;
    while (1) {
      if (*t2==0 || rc==0 ) break;
      rc *= (toupper(*t1++) == *t2++ ? 1 : 0);
    }
    if (rc==1) return (_values[jj]);
  }
  return NULL;
}

int StrBuffer::Separate(FILE *out) {
  const char mylim[]=" ";
  int sz = 0, cntr, num_fields;
  char *ptr;
  char *tmp = _buffer;
  
  while ( 1 ) { if (*tmp==0) break; if (*tmp==' ') sz++; tmp++; }
  if ( (sz-1)%2 != 0 ) {
    if ( out ) {
      fprintf(out, "Bummer: invalid variable definition. Should be \"variable = value\"\n");
    }
    return -1;
  }
  
  // check to see if we have enough memory, and grow the memory if needed
  sz = (sz-1)/2;
  _check_pair(sz);
  
  tmp = _buffer;
  cntr = 0; num_fields=0;
  ptr = strtok(tmp, mylim); 
  while ( ptr != NULL) {
    if ( num_fields == 0 ) {
      strcpy(_header, ptr);
    } else if ( (num_fields-1) % 2 == 0) {
      // put into keys
      strcpy(_keys[cntr], ptr);
    } else if ( (num_fields-1) % 2 != 0) {
      // put into values, and increase cntr
      strcpy(_values[cntr++], ptr);
    }
    num_fields++;
    ptr = strtok(NULL, mylim);
  }
  
  _pair_used = cntr;
  return cntr;
}

// Local Variables:
// mode: c++
// End:
