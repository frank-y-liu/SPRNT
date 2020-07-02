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

/**************************************************
  interface to the linear solvers
***************************************************/

#ifndef _SOLVER_INTERFACE_H
#define _SOLVER_INTERFACE_H

#include "solver_def.h"

class Solver_Interface {
 protected:
  static const char *_name;
  void        *_handle;
  void        *_s;

  create_solver_t create_solver;
  delete_solver_t delete_solver;
  init_solver_t   init_solver;
  solve_t         solve;
  get_type_t      get_type;
  clear_t         clear_solver;

public:
  Solver_Interface():_handle(0),_s(0) {}
  ~Solver_Interface();

  int Setup(int);
  void Init(int nnz, int dim, int *rows, int *cols, double *vals);
  int  Solve(double *rhs, double *x);
  int  SolverType();
  void Clear();
};

#endif

// Local Variables:
// mode: c++
// End:
