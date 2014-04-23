/*********************************************************************
  Copyright (C) 2011, 2014 International Business Machines
  All Rights Reserved

  Author:  Frank Liu, IBM
*********************************************************************/

/***************************************************** 
  interface for the dynamic library
******************************************************/

#include "lsolver.h"

extern "C" void* create_solver(int tp) {
  if (tp==L_UMF) {
    UMF* U=new UMF();
    return (void*)U;
  } else if ( tp==L_KLU) {
    KLU* K=new KLU();
    return (void*)K;
  } else {
    return ( NULL );
  }
}

extern "C" void delete_solver(void *s) {
  Solver *LS=(Solver*)s;
  if (LS->Type()==L_UMF) {
    UMF *U=(UMF*)s;
    delete U;
  } else {
    KLU *K=(KLU*)s;
    delete K;
  }
}

extern "C" void init_solver(void* s,int nnz, int dim, int* rows, int* cols, double* vals) {
  Solver *LS=(Solver*)s;
  LS->Init(nnz, dim, rows, cols, vals);
}

extern "C" void factor_solver(void *s) {
  Solver *LS=(Solver*)s;
  if (LS->Type() == L_UMF) {
    UMF* U=(UMF*)s;
    U->Factor();
  } else {
    KLU* K=(KLU*)s;
    K->Factor();
  }
}

extern "C" int solve(void *s, double *rhs, double *x) {
  Solver *LS=(Solver*)s;
  if (LS->Type() == L_UMF) {
    UMF *U=(UMF*)s;
    return (U->Solve(rhs, x));
  } else {
    KLU *K=(KLU*)s;
    return (K->Solve(rhs, x));
  }
}

extern "C" int get_type(void *s) {
  Solver *LS=(Solver*)s;
  return (LS->Type());
}

extern "C" void clear(void *s) {
  Solver *LS=(Solver*)s;
  LS->Clear();
}

// Local Variables:
// mode: c++
// End:
