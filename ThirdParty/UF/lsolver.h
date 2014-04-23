/*********************************************************************
  Copyright (C) 2011, 2014 International Business Machines
  All Rights Reserved

  Author:  Frank Liu, IBM
*********************************************************************/

/*************************************************** 
   interface to the linear solvers 
***************************************************/

#ifndef _LSOLVER_H
#define _LSOLVER_H

#include "solver_def.h"
#include "solver_vec.h"

#include "umfpack.h"
#include "klu.h"

/* Base solver class */
class	Solver {
 protected:

  int      _type;
  int      _is_factored;
  int      _nnz;
  int      _dim;

  int     *_rows_ptr;  // these are pointers, set up at init time
  int     *_cols_ptr;
  double  *_vals_ptr;

  flexvec<int>   _rows;
  flexvec<int>   _cols;
  flexvec<double> _vals;

 public:
  Solver():_type(L_UNKNOWN),_is_factored(0),_nnz(0),_dim(0),
	   _rows_ptr(NULL),_cols_ptr(NULL),_vals_ptr(NULL) {}
  virtual ~Solver() {}

  virtual void Init( int nnz, int dim, int *rows, int *cols, double *vals);
  virtual int Factor() { return (-1); }
  virtual int Solve(double *rhs, double *x) { return (-1); }
  virtual int Type() const { return _type; }
  virtual void Clear() { _is_factored=0; }
};

/* UMF: yet another general purpose direct solver */
class	UMF : public Solver  {
protected:
  int             _NRC;

  flexvec<int>    _Ap;
  flexvec<int>    _Ai;
  flexvec<double> _Ax;

  void       *_Symbolic;
  void       *_Numeric;

public:
  UMF() : Solver() {
    _type=L_UMF;_NRC=0;_Symbolic=NULL; _Numeric=NULL;
  }
  ~UMF() {
    if (_Symbolic) umfpack_di_free_symbolic(&_Symbolic);
    if (_Numeric)  umfpack_di_free_numeric( &_Numeric);
  }

  int Factor();
  int Solve(double *rhs, double *x);
};

/* KLU: general purpose solver for circuits */
class	KLU : public Solver  {
protected:
  int         _NRC;

  flexvec<int>     _Ap;
  flexvec<int>     _Ai;
  flexvec<double>  _Ax;

  void       *_Symbolic;
  void       *_Numeric;

public:
  KLU() : Solver() {
    _type=L_KLU;_NRC=0;_Symbolic=NULL; _Numeric=NULL;
  }
  ~KLU();

  int Factor();
  int Solve(double *rhs, double *x);

};


#endif

// Local Variables:
// mode: c++
// End:
