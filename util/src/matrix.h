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
  History
  09/10/11    FYL   Created 
*************************************************************************/

#ifndef _MATRIX
#define _MATRIX

#include <stdio.h>
#include <string.h>

#include "base.h"
#include "blaspr.h"
#include "flexvec.h"

/*
  Debugging facilities
 */
void matrix_dump(FILE *fp, char* mnm, int nr, int nc, double* h, int ldh );
void matrix_dump(FILE *fp, char* mnm, int nr, int nc, dComplex* h, int ldh );
void vector_dump(FILE *fp, char* mnm, int nr, double* v );
void vector_dump(FILE *fp, char* mnm, int nr, dComplex* v);


enum TransType { NOT_TRANS=0, TRANS, HERMIT};  // as the name says
enum WhereType { DOWN=0, RIGHT, DIAG };        // for concatenation 


inline char transchar( TransType trans ) {
  char   trc[] = "NTC";
  return trc[trans];
}

/* 
   Base class of matrix, defines the interface only
 */
class Matrix {
 protected:
  int     _dim[2];
 public: 
  Matrix() { _dim[0] = -1; }
  ~Matrix() {}
  int  Ndim(int i) const { return _dim[i]; };
  int  IsSized() const { return _dim[0] != -1; }; // dim[0] == -1 flags "uninitialized"
  void SetDimensions(int row,  int col)  { _dim[0]=row; _dim[1]=col; };
  void GetDimensions(int &row, int &col) const {  row=_dim[0]; col=_dim[1]; };
  void Clear() { _dim[0] = -1; }
};

/*
  Diagonal sparse matrix 
  
  Caveat: Only diagonal entries are stored. Any entrie not being "created" 
          will have garbage value. To avoid this side effect, consider "Zero" 
          all entries before using it (with some runtime penalties). Zero'ing
          can be avoided only if the user is certain all entries will be created. 

  Observation: a diagonal matrix does not have to be square

 */
class DiagMatrix : public Matrix {
  FlexVec<double> _mat;
 public:
  DiagMatrix() : Matrix(){}
  ~DiagMatrix() {}
  void CreateEntry(int ix, double val) { 
    _mat.grow(ix+1); _mat[ix] = val; _dim[0]=(_dim[0]>ix)?_dim[0]:ix+1; 
  }
  void SimplePrint(FILE *of);
  void MatlabDump(FILE *of, const char* Name="NoName");
  int Nnz() const { return MIN(_dim[0], _dim[1]); };
  double *Data()      { return &(_mat[0]);  };
  double &Data(int i) { return _mat[i]; };
  double &operator()(int i, int j) { assert(i==j); return _mat[i]; };
  void Zero(int rows, int cols) { 
    assert ( _dim[0] != -1 && rows >= 0 && cols >= 0);
    _dim[0] = rows; 
    _dim[1] = cols;
    _mat.size(Nnz());
    memset( _mat, 0, Nnz()*sizeof(double) ); 
  }
  void Eye(int rows, int cols, double alpha=1.0) { 
    assert ( _dim[0] != -1 && rows >= 0 && cols >= 0 );
    _dim[0] = rows;
    _dim[1] = cols;
    int n = Nnz(); const int zeroi=0; const int onei=1;
    _mat.size(n);
    FORTRAN(dcopy)( &n, &alpha, &zeroi, _mat, &onei );
  }
};

/*
  special selection matrix
 */
class SelMatrix : public Matrix {
  int          _sz;
  FlexVec<int> _sm;
public:
  SelMatrix() : Matrix(), _sz(0) {};
  void Clear() { Matrix::Clear(); _sz=0; }
  void AddEntry( int nd, int n ) { 
    _sm.grow( nd+1 ); 
    while (nd>=_sz) _sm[_sz++]=-1; 
    _sm[nd] = n;
  };
  void MatlabDump( FILE* fn, const char* nm );
  int  Select( int nd ) { if (nd<_sz) return _sm[nd]; else return(-1); } 
  int  Which( int n ) { 
    for( int nd=0; nd<_sz; ++nd )  
      if (_sm[nd] == n) return nd; 
    return -1; 
  } 
};

/*
  special adjacency matrix
 */
class AdjMatrix : public Matrix {
  int          _nr, _nc;
  FlexVec<int> _am;
public:
  AdjMatrix() : Matrix(), _nr(0), _nc(0) {};
  void Clear() { Matrix::Clear(); _nr=0; _nc=0; }
  void SetDimensions(int row,  int col)  { 
    assert(_nr == row && _nc <= col); 
    Matrix::SetDimensions(row, col); 
  } 
  int  Np(int k) { return _am[2*k]; }
  int  Nnz()     { return 2*Ndim(0); }
  int  Nn(int k) { return _am[2*k+1]; }
  int  AddBranch( int np, int nn ) { 
    _am.grow(2*_nr+2);
    _nc = MAX(np,nn,_nc);
    _am[2*_nr] = np; _am[2*_nr+1] = nn;
    return _nr++;
  }
  int  ForceBranch( int bn, int np, int nn ) {
    _nr=MAX(bn+1,_nr);
    _am.grow(2*_nr+2);
    _nc = MAX(np,nn,_nc);
    _am[2*bn] = np; _am[2*bn+1] = nn;
    return _nr;
  } 
  void MatlabDump( FILE* fn, const char* nm );
  void SimplePrint( FILE* fn );
};

/*
  General purpose sparse matrix
  
 */
class SparseMatrix : public Matrix {
 protected:
  FlexVec<int>    _irow;
  FlexVec<int>    _icol;
  FlexVec<double> _mat;
  FlexVec<int>    _crow;
  int     _is_sorted;  // flag
  int     _size;
  int     _asz;
 public:
  void MakeSpace(int ns) { 
    if (_asz < ns) { 
      _asz = _irow.grow(ns); 
      _asz = _icol.grow(ns); 
      _asz = _mat.grow(ns); 
    } 
  }
  void CreateEntryNC(int ix, int jx, double val) { 
    _irow[_size] = ix; 
    _icol[_size] = jx; 
    _mat[_size++] = val;
  };
  SparseMatrix() : Matrix(), _size(0) { _asz = _irow.size(0); } // set "uninitialized"
 ~SparseMatrix() {};
  int     Nnz() const { return _size; };
  int    *Row()  { return &(_irow[0]); };
  int    *Crow() { return &(_crow[0]); };
  int    *Col()  { return &(_icol[0]); }; 
  double *Data() { return &(_mat[0]);  };
  int    &Row(int i)  { return _irow[i]; };
  int    &Crow(int i) { return _crow[i]; };
  int    &Col(int i)  { return _icol[i]; };
  double &Data(int i) { return _mat[i]; };
  int  IsSetup() const { return IsSized() && _is_sorted; };
  void Clear() {_size = 0; _dim[0] = -1; };
  void CreateEntry(int ix, int jx, double val) { 
    MakeSpace(_size+1); _irow[_size] = ix; _icol[_size] = jx; _mat[_size] = val; ++_size; _is_sorted = 0; 
  };

  void Scale(double alpha) {
    const int onei = 1;
    assert( IsSized() );
    FORTRAN(dscal)(&_size, &alpha, Data(), &onei);
  }
  void Eye(int rows, int cols, double alpha=1.0);
  void Zero(int rows, int cols) { _size = 0; _is_sorted=0; _dim[0]=rows; _dim[1]=cols; };

  void SimplePrint(FILE *of);
  void MatlabDump(FILE *of, const char* Name="NoName");
};


#endif // _MATRIX

// Local Variables:
// mode: c++
// End:
