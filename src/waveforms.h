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
  simple class to store time series
 */

#ifndef _WAVEFORMS_H
#define _WAVEFORMS_H

#include "tdarray.h"
#include "flexvec.h"

/* 
   waveform class to hold the unsteady results
   built upon tdarray

   storage arrangement:
     height - nodes 
     width  - variables (Q, A, etc)
     length - time pts
     
   An example of usage:
    WV.Init( <number-of-nodes> ); // init

    rc = WV.Advance();           // advance before storing, 
    if (rc < 0) { something is wrong, go doodoo }
    WV.StoreQplusA(time, X);
    WV.StoreD(time, D_array);
    WV.StoreZ(time, Z_array);

    repeat
    
    use any of the query methods to get the results

    The same container can be reused by simply do 
    WV.Reset();

    each container can also be "cloned" by 
    Waveform N_WV(WV);

    but watch out the memory usage. If the number of nodes is large
    it will grow quite rapidly
 */

#define DFT_STORE_LEN 64 // enough for 1 hour at about 1 min interval

class Waveforms {
protected:
  enum      {Q_STORE=0, A_STORE, D_STORE, Z_STORE, UNK_STORE};

  int            _num_banks;   // number of banks, hard coded to 4
  int            _height;      // number of nodes
  int            _time_idx;    // current time point index
  
  ThreeDArray<double>                _data;  // data
  GrowVec<double,DFT_STORE_LEN>      _time;  // time points


public:
  Waveforms():_num_banks(0),_height(0),_time_idx(-1) {}
  Waveforms(int num_nodes):_num_banks(UNK_STORE),_height(num_nodes),
  			   _time_idx(-1),_data(UNK_STORE, num_nodes,DFT_STORE_LEN) {
    _time.size(DFT_STORE_LEN);
  }

  // copy constructor, copy everything!
  Waveforms( const Waveforms& other ):_num_banks(other._num_banks),
				      _height(other._height),
				      _time_idx(other._time_idx),
				      _data(other._data) {
    // except growvec doesn't have a copy constructor
    _time.size(_time_idx+1);
    for (int jj=0; jj<=_time_idx;jj++) _time[jj]=other._time.getVal(jj);
  }

  ~Waveforms() {}

  /* methods */

  // helper function to return the CLOSEST time index point
  // for a given value
  int FindTimeIdx(double t) {
    if (_time_idx < 0 ) return (-1);
    if (t < _time[0] ) return 0;
    if (t > _time[_time_idx]) return _time_idx;

    int jj, pf, cf;
    pf = ( t > _time[0] );
    // linear search
    for (jj=1; jj<=_time_idx; jj++) {
      cf = ( t > _time[jj] );
      if ( pf && !cf ) break;
      pf = cf;
    } 
    // has to be betwene jj-1 and jj
    if ( t-_time[jj-1] < _time[jj]-t ) return (jj-1);
    else return (jj);
  }

  void Init(int num_nodes) {
    _num_banks=UNK_STORE; // the last one in the entry
    _height = num_nodes;
    _time_idx = -1;
    _data.Init(UNK_STORE, num_nodes, DFT_STORE_LEN);
    _time.size(DFT_STORE_LEN);
  }

  // advance time point returns -1 when error
  int Advance() { _time_idx++; return (_data.AdvanceFrame()); }
  void Reset()  { _data.Reset(); _time_idx = -1; }

  // D and Z are derived
  void StoreD(double t, double *RHS) {
    _time[_time_idx]=t;
    _data.Store(_height, D_STORE, RHS); 
  }
  void StoreZ(double t, double *RHS) {
    _time[_time_idx] = t;
    _data.Store(_height, Z_STORE, RHS);
  }
  // Q and A can be directly extracted from X
  void StoreQplusA(double t, double *X) {
    const int onei=1;
    const int twoi=2;
    const int zeroi=0;

    _time[_time_idx] = t;
    // here is the black magic:
    //   Q are stored at 0, 2, 4 ..., A at 1, 3, 5
    _data.Store(_height, Q_STORE, zeroi, twoi, X);
    _data.Store(_height, A_STORE, onei,  twoi, X);
  }
  
  // access methods
  int Length() const { return (_time_idx+1); }
  int Width() const { return _num_banks; }
  int Height() const { return _height; }

  double FirstTimePt() const { return _time.getVal(0); }
  double LastTimePt() const { return _time.getVal(_time_idx); }
  double* TimeSeq() { return (double*)_time; }
  int     TimeSeq(double *target) { 
    for(int jj=0; jj<=_time_idx+1; jj++) target[jj]=_time.getVal(jj);
    return (_time_idx+1);
  }

  /* these methods could be more useful, 
     Note:  caller owns memory 
  */

  // get time varying Q, should be used in junction with TimeSeq()
  int GetQatNode(int n_idx, double *target) {
    assert ( n_idx < _height );
    return (_data.GetSlice(Q_STORE, n_idx, target)); 
  }

  int GetAatNode(int n_idx, double *target) {
    assert ( n_idx < _height );
    return (_data.GetSlice(A_STORE, n_idx, target)); 
  }

  int GetDatNode(int n_idx, double *target) {
    assert ( n_idx < _height );
    return (_data.GetSlice(D_STORE, n_idx, target)); 
  }

  int GetZatNode(int n_idx, double *target) {
    assert ( n_idx < _height );
    return (_data.GetSlice(Z_STORE, n_idx, target)); 
  }

  // get Q across all nodes at a given time point
  int GetQatTime(int t_idx, double *target) {
    assert ( t_idx <= _time_idx );
    return (_data.GetCol(t_idx, Q_STORE, target));
  }
  int GetAatTime(int t_idx, double *target) {
    assert ( t_idx <= _time_idx );
    return (_data.GetCol(t_idx, A_STORE, target));
  }
  int GetDatTime(int t_idx, double *target) {
    assert ( t_idx <= _time_idx );
    return (_data.GetCol(t_idx, D_STORE, target));
  }
  int GetZatTime(int t_idx, double *target) {
    assert ( t_idx <= _time_idx );
    return (_data.GetCol(t_idx, Z_STORE, target));
  }

  // get Q for a perticular node, at a time point
  double GetNodeQatTime(int t_idx, int n_idx) {
    assert ( t_idx <= _time_idx && n_idx < _height );
    double tmp[16];  // hmmm, a little hand-waivering
    _data.GetRow(t_idx,n_idx, &tmp[0]); 
    return (tmp[Q_STORE]);
  }

  double GetNodeAatTime(int t_idx, int n_idx) {
    assert ( t_idx <= _time_idx && n_idx < _height );
    double tmp[16]; 
    _data.GetRow(t_idx,n_idx, &tmp[0]); 
    return (tmp[A_STORE]);
  }

  double GetNodeDatTime(int t_idx, int n_idx) {
    assert ( t_idx <= _time_idx && n_idx < _height );
    double tmp[16]; 
    _data.GetRow(t_idx,n_idx, &tmp[0]); 
    return (tmp[D_STORE]);
  }

  double GetNodeZatTime(int t_idx, int n_idx) {
    assert ( t_idx <= _time_idx && n_idx < _height );
    double tmp[16]; 
    _data.GetRow(t_idx,n_idx, &tmp[0]); 
    return (tmp[Z_STORE]);
  }
  
};

#endif


// Local variables:
// mode: c++
// End:
