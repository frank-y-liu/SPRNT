/*********************************************************************
  Copyright (C) 2011, 2014 International Business Machines
  All Rights Reserved

  Author:  Frank Liu, IBM
*********************************************************************/

/***********************************************
 methods for lsolver class 
***********************************************/
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "lsolver.h"

void Solver::Init(int nnz, int dim, int *rows, int* cols, double *vals) {
  assert (nnz>0 && dim>0);
  assert( rows && cols && vals);

  _nnz = nnz;
  _dim = dim;

  _rows_ptr = rows;
  _cols_ptr = cols;
  _vals_ptr = vals;

  _rows.size(_nnz);
  _cols.size(_nnz);
  _vals.size(_nnz);

}

int UMF::Factor() {
  int status;
  double Info[UMFPACK_INFO];  // UMF internal communcation array
  double Control[UMFPACK_CONTROL];

  memcpy(&_rows[0], _rows_ptr, _nnz*sizeof(int));
  memcpy(&_cols[0], _cols_ptr, _nnz*sizeof(int));
  memcpy(&_vals[0], _vals_ptr, _nnz*sizeof(double));
  
  _Ap.size(_dim+1);
  _Ai.size(_nnz);
  _Ax.size(_nnz);

  _NRC = _nnz;

  /* get the default control parameters */
  umfpack_di_defaults (Control) ;
  
  /* convert to column form */
  status = umfpack_di_triplet_to_col(_dim, _dim, _NRC, 
				     &_rows[0], &_cols[0], 
				     &_vals[0], &_Ap[0], &_Ai[0], &_Ax[0], (int*)NULL) ;
  if (status < 0 ) {
    umfpack_di_report_status (Control, status) ;
    fprintf(stderr,"umfpack_di_triplet_to_col failed\n");
    return (-1);
  }
  
  /* ---------------------------------------------------------------------- */
  /* symbolic factorization */
  /* ---------------------------------------------------------------------- */
  
  if ( _Symbolic == NULL ) { // skip the symbolic factorization steps except the very first one
    status = umfpack_di_symbolic (_dim, _dim, &_Ap[0], &_Ai[0], &_Ax[0], &_Symbolic, Control, Info) ;
    if (status < 0 ) {
      umfpack_di_report_info (Control, Info) ;
      umfpack_di_report_status (Control, status) ;
      fprintf(stderr,"umfpack_di_symbolic failed\n");
      return (-1);
    }
  }
  
  /* ---------------------------------------------------------------------- */
  /* numeric factorization */
  /* ---------------------------------------------------------------------- */
  if ( _Numeric )  umfpack_di_free_numeric (&_Numeric);
  status = umfpack_di_numeric (&_Ap[0], &_Ai[0], &_Ax[0], _Symbolic, &_Numeric, Control, Info) ;
  if (status < 0 ) {
    umfpack_di_report_info (Control, Info) ;
    umfpack_di_report_status (Control, status) ;
    fprintf(stderr,"umfpack_di_numeric failed\n") ;
    return (-1);
  }
  _is_factored = 1;

  return 0;
}

int UMF::Solve(double *rhs, double *x) {
  int status;

  if (!_is_factored) status = Factor();
  
  int NEQN = _dim;

  double Info[UMFPACK_INFO], Control[UMFPACK_CONTROL];
  
  /* load the rhs only (matrix already there) */
  if (x != rhs) memcpy( x, rhs, NEQN*sizeof(double) ); 
  
  /* deal with transpose later */
  status = umfpack_di_solve (UMFPACK_A, &_Ap[0], &_Ai[0], &_Ax[0], x, rhs, _Numeric, Control, Info);
  
  return (status);
}

/***
  KLU 
 ***/
int KLU::Factor() {
  int status;
  int ndim;

  _NRC = _nnz;
  ndim = _dim;
  
  memcpy(&_rows[0], _rows_ptr, _nnz*sizeof(int));
  memcpy(&_cols[0], _cols_ptr, _nnz*sizeof(int));
  memcpy(&_vals[0], _vals_ptr, _nnz*sizeof(double));
  
  _Ap.size(_dim+1);
  _Ai.size(_nnz);
  _Ax.size(_nnz);

  status = umfpack_di_triplet_to_col (ndim,ndim, _NRC, &_rows[0], &_cols[0], &_vals[0], 
				      &_Ap[0], &_Ai[0], &_Ax[0], (int*)NULL) ;
  if (status < 0 ) {
    fprintf(stderr,"KLU: umfpack_di_triplet_to_col failed\n");
    return (-1);
  }

  klu_common Common;
  klu_defaults( &Common );
  klu_numeric *Num;

  if ( _Symbolic == NULL ) _Symbolic = (void*)klu_analyze(ndim, &_Ap[0], &_Ai[0], &Common);
  if ( _Symbolic == NULL ) {
    fprintf(stderr,"KLU: symbolic analysis failed\n");
    return (-1);
  }

  if ( _Numeric ) {
    Num = (klu_numeric*)_Numeric;
    klu_free_numeric( &Num, &Common);
  }

  _Numeric = (void*)klu_factor(&_Ap[0], &_Ai[0], &_Ax[0], (klu_symbolic*)_Symbolic, &Common);
  if ( _Numeric == NULL ) {
    fprintf(stderr,"KLU: numeric factorization failed\n");
    return (-1);
  }
  _is_factored = 1;

  return 0;
}

int KLU::Solve(double *rhs, double *x) {
  int rc = 0;

  if (!_is_factored) rc = Factor();
  
  int NEQN = _dim;
  
  klu_symbolic *Sym = (klu_symbolic*)_Symbolic;
  klu_numeric  *Num = (klu_numeric*)_Numeric;

  klu_common Common;
  klu_defaults( &Common );
  
  /* load the rhs only (matrix already there) */
  if (x != rhs) memcpy( x, rhs, NEQN*sizeof(double) ); 
  
  /* we will deal with transpose later */
  rc = klu_solve(Sym, Num, NEQN, 1, x, &Common);

  return rc;
}

KLU::~KLU() {
  klu_common Common;
  klu_defaults(&Common);
  klu_symbolic *Sym = (klu_symbolic*)_Symbolic;
  klu_numeric  *Num = (klu_numeric*)_Numeric;
  
  if ( _Symbolic ) klu_free_symbolic( &Sym, &Common);
  if ( _Numeric  ) klu_free_numeric ( &Num,  &Common);
}


// Local Variables:
// mode: c++
// End:
