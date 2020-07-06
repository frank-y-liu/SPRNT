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
 a map to match strings to integers. Insertion only, no deletions
 plus a few inquiry

 anticipated number of unique entries 100 ~ 100,000

 use unordered_map

 interface using char* to be backward compatible

 capacity: since we use int, the max number of entries is ~2,000,000,000

 not particularly fast, but we only use it when parsing the netlist.

 Update:
   on older compiler which doesn't support c++0x, use <map>
   unscientific testing shows <map> is twice as slow as <unordered_map>

 */

#ifndef _SMAP_H
#define _SMAP_H

#if __cplusplus > 201100L
#define HAS_CXX_SUPPORT
#endif

#include <string>
#ifdef HAS_CXX_SUPPORT
#include <unordered_map>
#else
#include <map>
#endif

using namespace std;

class SMap {
private:
#ifdef HAS_CXX_SUPPORT
  typedef std::unordered_map<std::string, int> MapType;
#else
  typedef std::map<std::string, int> MapType;
#endif
  int _cntr;
  MapType _M;

public:
#ifdef HAS_CXX_SUPPORT
  SMap() : _cntr(0) { _M.rehash(15000000); } // slightly different hash function
#else
  SMap() : _cntr(0) {}
#endif

  ~SMap() {}

  /* public methods */
  // return the size of the mapping table
  int Size() const { return (_M.size()); }

  // check, if the entry is found returns the key
  // if not, create a new entry and returns key
  int Create(const char *s) {
    string ns(s);
    MapType::const_iterator loc;

    loc = _M.find(ns);
    if (loc == _M.end()) {
      _M[ns] = _cntr;
      return (_cntr++);
    } else {
      return _M[ns];
    }
  }

  // check only
  // if not found, return -1
  // if found, return key
  int Check(const char *s) {
    string ns(s);
    MapType::const_iterator loc;

    loc = _M.find(ns);
    if (loc == _M.end()) {
      return (-1);
    } else {
      return _M[ns];
    }
  }
};

#undef HAS_CXX_SUPPORT
#endif

// Local Variables:
// mode: c++
// End:
