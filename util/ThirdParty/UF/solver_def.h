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

#ifndef _SOLVER_DEF_H
#define _SOLVER_DEF_H

#define L_UMF 0
#define L_KLU 1
#define L_UNKNOWN 2

typedef  void*(*create_solver_t)(int);
typedef  void (*delete_solver_t)(void*);
typedef  void (*init_solver_t)(void *, int, int, int*, int*, double*);
typedef  int  (*solve_t)(void*, double*, double*);
typedef  int  (*get_type_t)(void*);
typedef  void (*clear_t)(void*);

#endif

// Local Variables:
// mode: c++
// End:
