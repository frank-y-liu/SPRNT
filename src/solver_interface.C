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

/*************************************
  methods for the solver interface 
*******************************************/

#include <stdio.h>
#include <dlfcn.h>

#include "solver_interface.h"

const char* Solver_Interface::_name="libsolvers.so";

Solver_Interface::~Solver_Interface() {
  if ( _s ) delete_solver(_s);
  if (_handle) dlclose(_handle);
}

int Solver_Interface::Setup(int tp) {
  _handle = dlopen(_name, RTLD_LAZY);
  if (!_handle) {
    fprintf(stderr, "%s\n", dlerror());
    return (-1);
  }
  
  // function names and types hardcoded
  create_solver = (create_solver_t) dlsym(_handle,"create_solver");
  delete_solver = (delete_solver_t) dlsym(_handle, "delete_solver");
  init_solver = (init_solver_t) dlsym(_handle,"init_solver");
  solve = (solve_t)dlsym(_handle,"solve");
  get_type = (get_type_t)dlsym(_handle, "get_type");
  clear_solver = (clear_t)dlsym(_handle, "clear");

  // check to make sure all functions are available
  if (!create_solver) {
    fprintf(stderr,"%s\n", dlerror());
    return (-1);
  }
  if (!delete_solver) {
    fprintf(stderr,"%s\n", dlerror());
    return (-1);
  }
  if (!init_solver) {
    fprintf(stderr,"%s\n", dlerror());
    return (-1);
  }
  if (!solve) {
    fprintf(stderr,"%s\n", dlerror());
    return (-1);
  }
  if (!get_type) {
    fprintf(stderr,"%s\n", dlerror());
    return (-1);
  }
  if (!clear_solver) {
    fprintf(stderr,"%s\n", dlerror());
    return (-1);
  }

  // make sure we can allocate a solver
  _s = create_solver(tp);
  if (!_s) {
    fprintf(stderr,"Bummer, unable to allocate solver of type %d\n",tp);
    return (-1);
  }
  return 0;
}

void Solver_Interface::Init(int nnz, int dim, int *rows, int *cols, double *vals) {
  init_solver(_s, nnz, dim, rows, cols,vals);
}

int Solver_Interface::Solve(double *rhs, double *x) {
  return ( solve(_s, rhs, x) );
}

int Solver_Interface::SolverType() {
  return ( get_type(_s));
}

void Solver_Interface::Clear() {
  clear_solver(_s);
}

// Local Variables:
// mode: c++
// End:
