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
  methods for node class
  08/20/11    methods for nodes and equations
***********************************************************************/

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "node.h"
#include "options.h"

#define USEDIFF 0

/* we need some place to hold the initializer for the class */
int Node::_debug = 0;
int ArSrcEqn::_debug = 0;
int QrSrcEqn::_debug = 0;
int QrDepEqn::_debug = 0;
int MasStvEqn::_debug = 0;
double MasStvEqn::_spin_coef = 1.0;

// Methods related to Node
void Node::AddLocationInfo(double x, double y) {
  _x = x;
  _y = y;
}

void Node::AssignEntries(int q_idx, int a_idx) {
  assert(a_idx >= 0 && q_idx >= 0);

  _q_idx = q_idx;
  _a_idx = a_idx;
}

void Node::SetInitValues(double *X) {
  X[_q_idx] = _q0;
  X[_a_idx] = _a0;
}

/* kinematic routing using bisection method, returns A */
/* if s_r <= 0, returns -1.0
   upper level code should take care of it */
double Node::KinematicEstimate(double Q) {
  const double tol = 1e-6;
  const int num_iter = 50; // all hard coded
  const double p1 = 2 / 3.0;
  int jj;
  double coef = sqrt(_sr / _nsq); // = sqrt(s0)/n
  double a0, a1, at;
  double r, drda;
  double v0, v1, vt;

  if (_sr < 5e-6) { // sr could be negative or zero, make coef NaN, or zero
    return -1.0;
  }

  a0 = 1e-6;
  a1 = 0.01;
  const int slimit =
      25; // 2^25 * 0.01 ~= 3.3e+5, hopefully large enough for upper bound
  int found_ub = 0;
  for (jj = 0; jj < slimit; jj++) {
    _xs->GetHydroRadius(a1, r, drda);
    v1 = coef * pow(r, p1) * a1 - Q;
    if (v1 > 0) {
      found_ub = 1;
      break;
    }
    a1 *= 2.0;
  }
  if (found_ub == 0) {
    return -1.0; // failed to find ub, let up-level handle the fault
  }

  // we use the bisection method
  _xs->GetHydroRadius(a0, r, drda);
  v0 = coef * pow(r, p1) * a0 - Q;

  for (jj = 0; jj < num_iter; jj++) {
    if (ABS(a0 - a1) < tol) {
      at = v0;
      break;
    }
    at = 0.5 * (a1 + a0);
    _xs->GetHydroRadius(at, r, drda);
    vt = coef * pow(r, p1) * at - Q;

    if (ABS(vt) < tol)
      break;
    else if (vt > 0) {
      v1 = vt;
      a1 = at;
    } else {
      v0 = vt;
      a0 = at;
    }
  }

  // we are not checking number of iters
  // just returns the average of the last two iterations
  // results could be crappy if predefined iteration limit is not sufficient
  return (a0 + a1) / 2.0;
}

// calculate the froude number at a perticular node
// which is Q*sqrt(b)/ (sqrt_g* A^1.5)
//
// There is also a vectorized implementation of this function in the
// subcatchment class
double Node::GetFroude(double *X) {
  double q = X[QIdx()];
  double a = X[AIdx()];
  double b = this->_xs->GetWidth(a);

  assert(a > 0);
  assert(b > 0);

  return q * sqrt(b) / (OPT.SqrtG() * a * sqrt(a));
}

// overloaded with given q and a
double Node::GetFroude(double q, double a) {
  if (a < 0)
    return -1.0;

  double b = this->_xs->GetWidth(a);
  return (q * sqrt(b) / (OPT.SqrtG() * a * sqrt(a)));
}

/* Methods related to Equations and its derived classes */

/**************************************************************
  mass source term, the location is hard coded based on the type
***************************************************************/
void ArSrcEqn::CreateEntry(SparseMatrix *w, Node **nn) {
  w->CreateEntry(_row, nn[_nidx]->AIdx(), 0.0);
}

int ArSrcEqn::EvaluateRHS(Node **n, double *X, double *Xp, double *RHS) {
  int where = n[_nidx]->AIdx();
  double v = _S->Evaluate(0.0);

  RHS[_row] = X[where] - v;

  return OK;
}

int ArSrcEqn::Evaluate(Node **n, double *X, double *Xp, double *RHS,
                       SparseMatrix *M) {
  int where = n[_nidx]->AIdx();
  double v = _S->Evaluate(0.0); // this value should be the same as in the graph

  M->CreateEntryNC(_row, where, 1.0);
  RHS[_row] = X[where] - v;

  return OK;
}

int ArSrcEqn::EvaluateRHS(double t, double dt, Node **n, double *X, double *Xp,
                          double *RHS) {
  int where = n[_nidx]->AIdx();
  double v = _S->Evaluate(t);

  RHS[_row] = X[where] - v;

  return OK;
}

int ArSrcEqn::Evaluate(double t, double dt, Node **n, double *X, double *Xp,
                       double *RHS, SparseMatrix *M) {
  int where = n[_nidx]->AIdx();
  double v = _S->Evaluate(t);

  M->CreateEntryNC(_row, where, 1.0);
  RHS[_row] = X[where] - v;

  return OK;
}

int ArSrcEqn::EvaluateJac(double t, double dt, Node **n, double *X, double *Xp,
                          SparseMatrix *M) {
  // obviously not yet implemented. ditto for others
  return OK;
}

void ArSrcEqn::PrintValues(FILE *F, double t, double dt, Node **n, double *X,
                           double *Xp, double *RHS) {
  int where = n[_nidx]->AIdx();
  fprintf(F, "ASRC: (%4d %4d)->%6.2e || (%4d)->%+6.2e\n", _row, where, 1.0,
          _row, RHS[_row]);
}

/**************************************************************
      mass dependent term
***************************************************************/
int ArDepEqn::AddDep(int dep_idx, double ratio) {
  assert(dep_idx >= 0);

  _dep_idx = dep_idx;
  _ratio = ratio;

  return OK;
}

int ArDepEqn::EvaluateJac(double t, double dt, Node **n, double *X, double *Xp,
                          SparseMatrix *M) {
  return OK;
}

int ArDepEqn::EvaluateRHS(Node **n, double *X, double *Xp, double *RHS) {
  int idx1 = n[_nidx]->AIdx();
  int idx2 = n[_dep_idx]->AIdx();

  RHS[_row] = (X[idx1] - _ratio * X[idx2]) * OPT.KRatio();
  return OK;
}

int ArDepEqn::Evaluate(Node **n, double *X, double *Xp, double *RHS,
                       SparseMatrix *M) {
  int idx1 = n[_nidx]->AIdx();
  int idx2 = n[_dep_idx]->AIdx();

  M->CreateEntryNC(_row, idx1, 1.0);
  M->CreateEntryNC(_row, idx2, -_ratio);
  RHS[_row] = (X[idx1] - _ratio * X[idx2]) * OPT.KRatio();

  return OK;
}

int ArDepEqn::EvaluateRHS(double t, double dt, Node **n, double *X, double *Xp,
                          double *RHS) {
  int idx1 = n[_nidx]->AIdx();
  int idx2 = n[_dep_idx]->AIdx();

  RHS[_row] = (X[idx1] - _ratio * X[idx2]) * OPT.KRatio();
  return OK;
}

int ArDepEqn::Evaluate(double t, double dt, Node **n, double *X, double *Xp,
                       double *RHS, SparseMatrix *M) {
  int idx1 = n[_nidx]->AIdx();
  int idx2 = n[_dep_idx]->AIdx();

  M->CreateEntryNC(_row, idx1, 1.0);
  M->CreateEntryNC(_row, idx2, -_ratio);
  RHS[_row] = (X[idx1] - _ratio * X[idx2]) * OPT.KRatio();

  return OK;
}

void ArDepEqn::PrintValues(FILE *F, double t, double dt, Node **n, double *X,
                           double *Xp, double *RHS) {
  int idx1 = n[_nidx]->AIdx();
  int idx2 = n[_dep_idx]->AIdx();

  fprintf(F, "ADEP: (%4d %4d)->%6.2e (%4d %4d)->%6.2e || (4%d)->%+6.2e\n", _row,
          idx1, 1.0, _row, idx2, -_ratio, _row, RHS[_row]);
}

void ArDepEqn::CreateEntry(SparseMatrix *w, Node **nn) {
  w->CreateEntry(_row, nn[_nidx]->AIdx(), 0.0);
  w->CreateEntry(_row, nn[_dep_idx]->AIdx(), 0.0);
}

/**************************************************************
 continuity equation in stv
  if _pp is presented, add the lateral flow term
***************************************************************/
int MasStvEqn::EvaluateRHS(Node **n, double *X, double *Xp, double *RHS) {
  // since we don't rely on continuity equation for Q propagation
  // this is degenerated. We simply create entries and make sure the matrix
  // is not singular
  RHS[_row] = 0.0;

  return OK;
}

int MasStvEqn::Evaluate(Node **n, double *X, double *Xp, double *RHS,
                        SparseMatrix *M) {
  int idx_up_a = n[_up_idx]->AIdx();
  int idx_up_q = n[_up_idx]->QIdx();
  int idx_dn_a = n[_dn_idx]->AIdx();
  int idx_dn_q = n[_dn_idx]->QIdx();

  RHS[_row] = 0;

  M->CreateEntryNC(_row, idx_up_a, 0);
  M->CreateEntryNC(_row, idx_dn_a, 0);
  M->CreateEntryNC(_row, idx_up_q, 1.0);
  M->CreateEntryNC(_row, idx_dn_q, 1.0);

  return OK;
}

int MasStvEqn::EvaluateRHS(double t, double dt, Node **n, double *X, double *Xp,
                           double *RHS) {
  double lambda2 = 2.0 * dt / _dx;
  double pv;
  int idx_up_a = n[_up_idx]->AIdx();
  int idx_up_q = n[_up_idx]->QIdx();
  int idx_dn_a = n[_dn_idx]->AIdx();
  int idx_dn_q = n[_dn_idx]->QIdx();

  RHS[_row] = X[idx_up_a] + X[idx_dn_a] - Xp[idx_up_a] - Xp[idx_dn_a] +
              lambda2 * (X[idx_dn_q] - X[idx_up_q]);

  if (_pp) {
    pv = _pp->Evaluate(t) * MasStvEqn::_spin_coef;
    RHS[_row] -= lambda2 * pv;
    if (OPT.DebugLevel() >= 2 && MasStvEqn::Debug() == 1) {
      printf("    Qstv %lx, PourPort(%.2e)=%.2e\n", (unsigned long)this, t, pv);
    }
  }
  return OK;
}

int MasStvEqn::Evaluate(double t, double dt, Node **n, double *X, double *Xp,
                        double *RHS, SparseMatrix *M) {
  double lambda2 = 2.0 * dt / _dx;
  double pv;
  int idx_up_a = n[_up_idx]->AIdx();
  int idx_up_q = n[_up_idx]->QIdx();
  int idx_dn_a = n[_dn_idx]->AIdx();
  int idx_dn_q = n[_dn_idx]->QIdx();

  RHS[_row] = X[idx_up_a] + X[idx_dn_a] - Xp[idx_up_a] - Xp[idx_dn_a] +
              lambda2 * (X[idx_dn_q] - X[idx_up_q]);

  if (_pp) {
    pv = _pp->Evaluate(t) * MasStvEqn::_spin_coef;
    RHS[_row] -= lambda2 * pv;
    if (OPT.DebugLevel() >= 2 && MasStvEqn::Debug() == 1) {
      printf("    Qstv %lx, PourPort(%.2e)=%.2e\n", (unsigned long)this, t, pv);
    }
  }

  M->CreateEntryNC(_row, idx_up_a, 1.0);
  M->CreateEntryNC(_row, idx_dn_a, 1.0);
  M->CreateEntryNC(_row, idx_up_q, -lambda2);
  M->CreateEntryNC(_row, idx_dn_q, lambda2);

  return OK;
}

int MasStvEqn::EvaluateJac(double t, double dt, Node **n, double *X, double *Xp,
                           SparseMatrix *M) {
  return OK;
}

void MasStvEqn::PrintValues(FILE *F, double t, double dt, Node **n, double *X,
                            double *Xp, double *RHS) {
  double lambda2 = 2.0 * dt / _dx;
  int idx_up_a = n[_up_idx]->AIdx();
  int idx_up_q = n[_up_idx]->QIdx();
  int idx_dn_a = n[_dn_idx]->AIdx();
  int idx_dn_q = n[_dn_idx]->QIdx();

  fprintf(F, "** Qup=%6.2e Aup=%6.2e Qdown=%6.2e Adown=%6.2e\n", X[idx_up_q],
          X[idx_up_a], X[idx_dn_q], X[idx_dn_a]);
  fprintf(F,
          "MASV: (%4d %4d)=%6.2e (%4d %4d)=%6.2e (%4d %4d)=%6.2e (%4d "
          "%4d)=%6.2e || (%4d)->%+6.2e\n",
          _row, idx_up_a, 1.0, _row, idx_dn_a, 1.0, _row, idx_up_q, -lambda2,
          _row, idx_dn_q, lambda2, _row, RHS[_row]);
}

void MasStvEqn::CreateEntry(SparseMatrix *w, Node **nn) {
  w->CreateEntry(_row, nn[_up_idx]->QIdx(), 0.0);
  w->CreateEntry(_row, nn[_up_idx]->AIdx(), 0.0);
  w->CreateEntry(_row, nn[_dn_idx]->QIdx(), 0.0);
  w->CreateEntry(_row, nn[_dn_idx]->AIdx(), 0.0);
}

/**************************************************************
   Q source term
***************************************************************/
int QrSrcEqn::EvaluateRHS(Node **n, double *X, double *Xp, double *RHS) {
  RHS[_row] = 0.0;
  return OK;
}

int QrSrcEqn::Evaluate(Node **n, double *X, double *Xp, double *RHS,
                       SparseMatrix *M) {
  int where = n[_nidx]->QIdx();
  M->CreateEntryNC(_row, where, 1.0);
  RHS[_row] = 0.0;

  return OK;
}

int QrSrcEqn::EvaluateRHS(double t, double dt, Node **n, double *X, double *Xp,
                          double *RHS) {
  int where = n[_nidx]->QIdx();
  double v = _S->Evaluate(t);

  RHS[_row] = X[where] - v;
  if (OPT.DebugLevel() >= 2 && QrSrcEqn::Debug() == 1) {
    printf("    Qsrc %lx, Y(%.2e)=%.2e\n", (unsigned long)this, t, v);
  }

  return OK;
}

int QrSrcEqn::Evaluate(double t, double dt, Node **n, double *X, double *Xp,
                       double *RHS, SparseMatrix *M) {
  int where = n[_nidx]->QIdx();
  double v = _S->Evaluate(t);

  M->CreateEntryNC(_row, where, 1.0);
  RHS[_row] = X[where] - v;
  if (OPT.DebugLevel() >= 2 && QrSrcEqn::Debug() == 1) {
    printf("    Qsrc %lx, Y(%.2e)=%.2e\n", (unsigned long)this, t, v);
  }

  return OK;
}

int QrSrcEqn::EvaluateJac(double t, double dt, Node **n, double *X, double *Xp,
                          SparseMatrix *M) {
  return OK;
}

void QrSrcEqn::PrintValues(FILE *F, double t, double dt, Node **n, double *X,
                           double *Xp, double *RHS) {
  int where = n[_nidx]->QIdx();

  fprintf(F, "QSRC: (%4d %4d)->%6.2e || (%4d)->%+6.2e\n", _row, where, 1.0,
          _row, RHS[_row]);
}

void QrSrcEqn::CreateEntry(SparseMatrix *w, Node **nn) {
  w->CreateEntry(_row, nn[_nidx]->QIdx(), 0.0);
}

/**************************************************************
   dependent Q terms
***************************************************************/
int QrDepEqn::AddDep(int sz, int *indices) {
  assert(sz > 0);

  _deps.size(sz);
  for (int j = 0; j < sz; j++)
    _deps[j] = indices[j];
  _num_deps = sz;
  return sz;
}

int QrDepEqn::EvaluateRHS(Node **n, double *X, double *Xp, double *RHS) {
  RHS[_row] = 0.0;
  return OK;
}

int QrDepEqn::Evaluate(Node **n, double *X, double *Xp, double *RHS,
                       SparseMatrix *M) {
  int idx1 = n[_nidx]->QIdx();
  int idx2;

  M->CreateEntryNC(_row, idx1, 1.0);
  for (int j = 0; j < _num_deps; j++) {
    idx2 = n[_deps[j]]->QIdx();
    M->CreateEntryNC(_row, idx2, -1.0);
  }
  RHS[_row] = 0.0;
  return OK;
}

/* if a _pp (pour port) is present, add the flow */
int QrDepEqn::EvaluateRHS(double t, double dt, Node **n, double *X, double *Xp,
                          double *RHS) {
  int idx1 = n[_nidx]->QIdx();
  int idx2;
  double v = X[idx1];
  double pv;

  for (int j = 0; j < _num_deps; j++) {
    idx2 = n[_deps[j]]->QIdx();
    v -= X[idx2];
  }
  RHS[_row] = v;
  if (_pp) {
    pv = _pp->Evaluate(t);
    RHS[_row] -= pv;
    if (OPT.DebugLevel() >= 2 && QrDepEqn::Debug() == 1) {
      printf("    Qdep %lx, Lateral Flow (%.2e)=%.2e\n", (unsigned long)this, t,
             v);
    }
  }
  /* we can "emphasize" the connections by scale it up in the rhs, it seems not
     necessary. Considering removing it */
  RHS[_row] *= OPT.KRatio();

  return OK;
}

int QrDepEqn::Evaluate(double t, double dt, Node **n, double *X, double *Xp,
                       double *RHS, SparseMatrix *M) {
  int idx1 = n[_nidx]->QIdx();
  int idx2;
  double v = X[idx1];
  double pv;

  M->CreateEntryNC(_row, idx1, 1.0);

  for (int j = 0; j < _num_deps; j++) {
    idx2 = n[_deps[j]]->QIdx();
    M->CreateEntryNC(_row, idx2, -1.0);
    v -= X[idx2];
  }
  RHS[_row] = v;
  if (_pp) {
    pv = _pp->Evaluate(t);
    RHS[_row] -= pv;
    if (OPT.DebugLevel() >= 2 && QrDepEqn::Debug() == 1) {
      printf("    Qdep %lx, PourPort(%.2e)=%.2e\n", (unsigned long)this, t, v);
    }
  }

  RHS[_row] *= OPT.KRatio();

  return OK;
}

int QrDepEqn::EvaluateJac(double t, double dt, Node **n, double *X, double *Xp,
                          SparseMatrix *M) {
  // to be implemented, see Evaluate() method
  return OK;
}

void QrDepEqn::PrintValues(FILE *F, double t, double dt, Node **n, double *X,
                           double *Xp, double *RHS) {
  int idx1 = n[_nidx]->QIdx();
  int idx2, j;

  fprintf(F, "QDEP: (%4d %4d)->%6.2e ", _row, idx1, 1.0);
  for (j = 0; j < _num_deps; j++) {
    idx2 = n[_deps[j]]->QIdx();
    fprintf(F, "(%4d %4d)->%6.2e ", _row, idx2, -1.0);
  }
  fprintf(F, "|| (%4d)->%+6.2e\n", _row, RHS[_row]);
}

void QrDepEqn::CreateEntry(SparseMatrix *w, Node **nn) {
  w->CreateEntry(_row, nn[_nidx]->QIdx(), 0.0);
  for (int j = 0; j < _num_deps; j++) {
    w->CreateEntry(_row, nn[_deps[j]]->QIdx(), 0.0);
  }
}

/**************************************************************
  momentum term in st venant equation
***************************************************************/
// The issue with this implementation is that each node is evaluated TWICE!!!
int MomStvEqn::Evaluate(Node **n, double *X, double *Xp, double *RHS,
                        SparseMatrix *M) {
  double half_g = OPT.Gravity() / 2.0;

  Node *UP = n[_up_idx];
  Node *DN = n[_dn_idx];

  int idx_up_a = UP->AIdx();
  int idx_up_q = UP->QIdx();
  int idx_dn_a = DN->AIdx();
  int idx_dn_q = DN->QIdx();

  double Qupsq = X[idx_up_q] * X[idx_up_q];
  double Qdnsq = X[idx_dn_q] * X[idx_dn_q];
  double Qupsqdiva = Qupsq / X[idx_up_a];
  double Qdnsqdiva = Qdnsq / X[idx_dn_a];

  double Qupsqabs = X[idx_up_q] * ABS(X[idx_up_q]);
  double Qdnsqabs = X[idx_dn_q] * ABS(X[idx_dn_q]);

  double Efdn = DN->XS()->GetEqFriction(X[idx_dn_a] + OPT.EpsilonA());
  double Efup = UP->XS()->GetEqFriction(X[idx_up_a] + OPT.EpsilonA());

  double up_n_corr = 1.0;
  double dn_n_corr = 1.0;
  double up_a = X[idx_up_a];
  double dn_a = X[idx_dn_a];

  if (OPT.HT() > 0 && up_a < 2.0 * OPT.HT()) {
    up_n_corr = OPT.HT() / up_a;
    if (up_n_corr < 10)
      up_n_corr = 1 + (log(1 + exp(12.0 * (up_n_corr - 1)))) / 12.0;
    if (up_n_corr > 20)
      up_n_corr = 20;
    if (OPT.DebugLevel() >= 2) {
      fprintf(stdout, "   HT invoked for a=%.2e, correction = %.2e\n", up_a,
              up_n_corr);
    }
  }

  if (OPT.HT() > 0 && dn_a < 2.0 * OPT.HT()) {
    dn_n_corr = OPT.HT() / dn_a;
    if (dn_n_corr < 10)
      dn_n_corr = 1 + (log(1 + exp(12.0 * (dn_n_corr - 1)))) / 12.0;
    if (dn_n_corr > 20)
      dn_n_corr = 20;
    if (OPT.DebugLevel() >= 2) {
      fprintf(stdout, "   HT invoked for a=%.2e, correction = %.2e\n", dn_a,
              dn_n_corr);
    }
  }

  RHS[_row] =
      (OPT.Beta() / _dx) * (Qdnsqdiva - Qupsqdiva) +
      (half_g / _dx) * (X[idx_dn_a] + X[idx_up_a]) *
          (DN->GetDepth(X[idx_dn_a]) - UP->GetDepth(X[idx_up_a])) -
      half_g * (DN->GetSR() * (X[idx_dn_a]) + UP->GetSR() * (X[idx_up_a])) +
      half_g * (DN->Nsq() * dn_n_corr * Qdnsqabs * Efdn +
                UP->Nsq() * up_n_corr * Qupsqabs * Efup);

  double up_fric_da = UP->XS()->GetEqFrictiondA(X[idx_up_a] + OPT.EpsilonA());
  double dn_fric_da = DN->XS()->GetEqFrictiondA(X[idx_dn_a] + OPT.EpsilonA());

  double dmdaup =
      (OPT.Beta() / _dx) * Qupsqdiva / (X[idx_up_a]) +
      (half_g / _dx) *
          (DN->GetDepth(X[idx_dn_a]) - UP->GetDepth(X[idx_up_a]) -
           (X[idx_up_a] + X[idx_dn_a]) * UP->GetDepthdA(X[idx_up_a])) -
      half_g * UP->GetSR() +
      half_g * UP->Nsq() * up_n_corr * Qupsqabs * up_fric_da;

  double dmdadn =
      (-OPT.Beta() / _dx) * Qdnsqdiva / (X[idx_dn_a]) +
      (half_g / _dx) *
          (DN->GetDepth(X[idx_dn_a]) - UP->GetDepth(X[idx_up_a]) +
           (X[idx_up_a] + X[idx_dn_a]) * DN->GetDepthdA(X[idx_dn_a])) -
      half_g * DN->GetSR() +
      half_g * DN->Nsq() * dn_n_corr * Qdnsqabs * dn_fric_da;

  M->CreateEntryNC(_row, UP->AIdx(), dmdaup);
  M->CreateEntryNC(_row, DN->AIdx(), dmdadn);
  M->CreateEntryNC(_row, UP->QIdx(), 0.0);
  M->CreateEntryNC(_row, DN->QIdx(), 0.0);

  return OK;
}

int MomStvEqn::EvaluateRHS(Node **n, double *X, double *Xp, double *RHS) {
  double half_g = OPT.Gravity() / 2.0;

  Node *UP = n[_up_idx];
  Node *DN = n[_dn_idx];

  int idx_up_a = UP->AIdx();
  int idx_up_q = UP->QIdx();
  int idx_dn_a = DN->AIdx();
  int idx_dn_q = DN->QIdx();

  double Qupsq = X[idx_up_q] * X[idx_up_q];
  double Qdnsq = X[idx_dn_q] * X[idx_dn_q];
  double Qupsqdiva = Qupsq / X[idx_up_a];
  double Qdnsqdiva = Qdnsq / X[idx_dn_a];

  double Qupsqabs = X[idx_up_q] * ABS(X[idx_up_q]);
  double Qdnsqabs = X[idx_dn_q] * ABS(X[idx_dn_q]);

  double Efdn = DN->XS()->GetEqFriction(X[idx_dn_a] + OPT.EpsilonA());
  double Efup = UP->XS()->GetEqFriction(X[idx_up_a] + OPT.EpsilonA());

  double up_n_corr = 1.0;
  double dn_n_corr = 1.0;
  double up_a = X[idx_up_a];
  double dn_a = X[idx_dn_a];

  if (OPT.HT() > 0 && up_a < 2.0 * OPT.HT()) {
    up_n_corr = OPT.HT() / up_a;
    if (up_n_corr < 10)
      up_n_corr = 1 + (log(1 + exp(12.0 * (up_n_corr - 1)))) / 12.0;
    if (up_n_corr > 20)
      up_n_corr = 20;
    if (OPT.DebugLevel() >= 2) {
      fprintf(stdout, "   HT invoked for a=%.2e, correction = %.2e\n", up_a,
              up_n_corr);
    }
  }

  if (OPT.HT() > 0 && dn_a < 2.0 * OPT.HT()) {
    dn_n_corr = OPT.HT() / dn_a;
    if (dn_n_corr < 10)
      dn_n_corr = 1 + (log(1 + exp(12.0 * (dn_n_corr - 1)))) / 12.0;
    if (dn_n_corr > 20)
      dn_n_corr = 20;
    if (OPT.DebugLevel() >= 2) {
      fprintf(stdout, "   HT invoked for a=%.2e, correction = %.2e\n", dn_a,
              dn_n_corr);
    }
  }

  RHS[_row] =
      (OPT.Beta() / _dx) * (Qdnsqdiva - Qupsqdiva) +
      (half_g / _dx) * (X[idx_dn_a] + X[idx_up_a]) *
          (DN->GetDepth(X[idx_dn_a]) - UP->GetDepth(X[idx_up_a])) -
      half_g * (DN->GetSR() * (X[idx_dn_a]) + UP->GetSR() * (X[idx_up_a])) +
      half_g * (DN->Nsq() * dn_n_corr * Qdnsqabs * Efdn +
                UP->Nsq() * up_n_corr * Qupsqabs * Efup);

  return OK;
}

int MomStvEqn::EvaluateJac(double t, double dt, Node **n, double *X, double *Xp,
                           SparseMatrix *M) {
  return OK;
}

int MomStvEqn::EvaluateRHS(double t, double dt, Node **n, double *X, double *Xp,
                           double *RHS) {
  double lambda2 = 2.0 * dt / _dx;
  double gdt = OPT.Gravity() * dt;

  Node *UP = n[_up_idx];
  Node *DN = n[_dn_idx];

  int idx_up_a = UP->AIdx();
  int idx_up_q = UP->QIdx();
  int idx_dn_a = DN->AIdx();
  int idx_dn_q = DN->QIdx();

  double Qupsq = X[idx_up_q] * X[idx_up_q];
  double Qdnsq = X[idx_dn_q] * X[idx_dn_q];

  double Qupsqdiva = Qupsq / X[idx_up_a];
  double Qdnsqdiva = Qdnsq / X[idx_dn_a];

  double Qupsqabs = X[idx_up_q] * ABS(X[idx_up_q]);
  double Qdnsqabs = X[idx_dn_q] * ABS(X[idx_dn_q]);

  double Efdn = DN->XS()->GetEqFriction(X[idx_dn_a] + OPT.EpsilonA());
  double Efup = UP->XS()->GetEqFriction(X[idx_up_a] + OPT.EpsilonA());

  double up_n_corr = 1.0;
  double dn_n_corr = 1.0;
  double up_a, dn_a;

  up_a = X[idx_up_a];
  dn_a = X[idx_dn_a];

  if (OPT.HT() > 0 && up_a < 2.0 * OPT.HT()) {
    up_n_corr = OPT.HT() / up_a;

    // better way, but more expensive
    if (up_n_corr < 10)
      up_n_corr = 1 + (log(1 + exp(12.0 * (up_n_corr - 1)))) / 12.0;
    if (up_n_corr > 20)
      up_n_corr = 20;
    if (OPT.DebugLevel() >= 2) {
      fprintf(stdout, "   HT invoked for a=%.2e, correction = %.2e\n", up_a,
              up_n_corr);
    }
  }

  if (OPT.HT() > 0 && dn_a < 2.0 * OPT.HT()) {
    dn_n_corr = OPT.HT() / dn_a;
    if (dn_n_corr < 10)
      dn_n_corr = 1 + (log(1 + exp(12.0 * (dn_n_corr - 1)))) / 12.0;
    if (dn_n_corr > 20)
      dn_n_corr = 20;
    if (OPT.DebugLevel() >= 2) {
      fprintf(stdout, "   HT invoked for a=%.2e, correction = %.2e\n", dn_a,
              dn_n_corr);
    }
  }

  // experimenting diffusive
  if (up_a < 1.5 * UP->XS()->GetMinA()) {
  }
#if USEDIFF
  double ff =
      1 / (1 + exp(-75.0 *
                   (X[idx_dn_a] - 0.1))); // hard coded, threshold=0.1, rate-75
  printf("a=%.4e ff=%.4e\n", X[idx_dn_a], ff);
  RHS[_row] =
      ff * (X[idx_dn_q] + X[idx_up_q] - Xp[idx_dn_q] - Xp[idx_up_q]) +
      lambda2 * (ff * OPT.Beta() * (Qdnsqdiva - Qupsqdiva) +
                 OPT.Gravity() * (X[idx_dn_a] + X[idx_up_a]) *
                     (DN->GetDepth(X[idx_dn_a]) - UP->GetDepth(X[idx_up_a]))) -
      gdt * (DN->GetSR() * (X[idx_dn_a]) + UP->GetSR() * (X[idx_up_a])) +
      gdt * (DN->Nsq() * dn_n_corr * Qdnsqabs * Efdn +
             UP->Nsq() * up_n_corr * Qupsqabs * Efup);
#else
  RHS[_row] =
      X[idx_dn_q] + X[idx_up_q] - Xp[idx_dn_q] - Xp[idx_up_q] +
      lambda2 * (OPT.Beta() * (Qdnsqdiva - Qupsqdiva) +
                 0.5 * OPT.Gravity() * (X[idx_dn_a] + X[idx_up_a]) *
                     (DN->GetDepth(X[idx_dn_a]) - UP->GetDepth(X[idx_up_a]))) -
      gdt * (DN->GetSR() * (X[idx_dn_a]) + UP->GetSR() * (X[idx_up_a])) +
      gdt * (DN->Nsq() * dn_n_corr * Qdnsqabs * Efdn +
             UP->Nsq() * up_n_corr * Qupsqabs * Efup);

#endif
  return OK;
}

// The issue with this implementation is that each node is evaluated TWICE!!!
int MomStvEqn::Evaluate(double t, double dt, Node **n, double *X, double *Xp,
                        double *RHS, SparseMatrix *M) {
  // this function can be CPU intensive, some optimization should be done!!!
  double lambda2 = 2.0 * dt / _dx;
  double gdt = OPT.Gravity() * dt;
  Node *UP = n[_up_idx];
  Node *DN = n[_dn_idx];

  int idx_up_a = UP->AIdx();
  int idx_up_q = UP->QIdx();
  int idx_dn_a = DN->AIdx();
  int idx_dn_q = DN->QIdx();

  double Qupsq = X[idx_up_q] * X[idx_up_q];
  double Qdnsq = X[idx_dn_q] * X[idx_dn_q];

  double Qupsqdiva = (Qupsq / X[idx_up_a]);
  double Qdnsqdiva = (Qdnsq / X[idx_dn_a]);

  double Qupsqabs = X[idx_up_q] * ABS(X[idx_up_q]);
  double Qdnsqabs = X[idx_dn_q] * ABS(X[idx_dn_q]);

  double Efdn = DN->XS()->GetEqFriction(X[idx_dn_a] + OPT.EpsilonA());
  double Efup = UP->XS()->GetEqFriction(X[idx_up_a] + OPT.EpsilonA());

  double up_n_corr = 1.0;
  double dn_n_corr = 1.0;
  double up_a, dn_a;

  up_a = X[idx_up_a];
  dn_a = X[idx_dn_a];

  if (OPT.HT() > 0 && up_a < 2.0 * OPT.HT()) {

    up_n_corr = OPT.HT() / up_a;
    if (up_n_corr < 10)
      up_n_corr = 1 + (log(1 + exp(12.0 * (up_n_corr - 1)))) / 12.0;
    if (up_n_corr > 20)
      up_n_corr = 20;
    if (OPT.DebugLevel() >= 2) {
      fprintf(stdout, "   HT invoked for a=%.2e, correction = %.2e\n", up_a,
              up_n_corr);
    }
  }

  if (OPT.HT() > 0 && dn_a < 2.0 * OPT.HT()) {
    dn_n_corr = OPT.HT() / dn_a;
    if (dn_n_corr < 10)
      dn_n_corr = 1 + (log(1 + exp(12.0 * (dn_n_corr - 1)))) / 12.0;
    if (dn_n_corr > 20)
      dn_n_corr = 20;
    if (OPT.DebugLevel() >= 2) {
      fprintf(stdout, "   HT invoked for a=%.2e, correction = %.2e\n", dn_a,
              dn_n_corr);
    }
  }
#if USEDIFF
  double ff =
      1 / (1 + exp(-75.0 *
                   (X[idx_dn_a] - 0.1))); // hard coded, threshold=0.1, rate-75
  RHS[_row] =
      ff * (X[idx_dn_q] + X[idx_up_q] - Xp[idx_dn_q] - Xp[idx_up_q]) +
      lambda2 * (ff * OPT.Beta() * (Qdnsqdiva - Qupsqdiva) +
                 OPT.Gravity() * (X[idx_dn_a] + X[idx_up_a]) *
                     (DN->GetDepth(X[idx_dn_a]) - UP->GetDepth(X[idx_up_a]))) -
      gdt * (DN->GetSR() * (X[idx_dn_a]) + UP->GetSR() * (X[idx_up_a])) +
      gdt * (DN->Nsq() * dn_n_corr * Qdnsqabs * Efdn +
             UP->Nsq() * up_n_corr * Qupsqabs * Efup);
#else
  RHS[_row] =
      X[idx_dn_q] + X[idx_up_q] - Xp[idx_dn_q] - Xp[idx_up_q] +
      lambda2 * (OPT.Beta() * (Qdnsqdiva - Qupsqdiva) +
                 0.5 * OPT.Gravity() * (X[idx_dn_a] + X[idx_up_a]) *
                     (DN->GetDepth(X[idx_dn_a]) - UP->GetDepth(X[idx_up_a]))) -
      gdt * (DN->GetSR() * (X[idx_dn_a]) + UP->GetSR() * (X[idx_up_a])) +
      gdt * (DN->Nsq() * dn_n_corr * Qdnsqabs * Efdn +
             UP->Nsq() * up_n_corr * Qupsqabs * Efup);
#endif

  double up_fric_da = UP->XS()->GetEqFrictiondA(X[idx_up_a] + OPT.EpsilonA());
  double dn_fric_da = DN->XS()->GetEqFrictiondA(X[idx_dn_a] + OPT.EpsilonA());

  double dmdqup = 1.0 -
                  2.0 * lambda2 * OPT.Beta() * (X[idx_up_q]) / (X[idx_up_a]) +
                  2.0 * gdt * UP->Nsq() * up_n_corr * Efup * ABS(X[idx_up_q]);
  double dmdqdn = 1.0 +
                  2.0 * lambda2 * OPT.Beta() * (X[idx_dn_q]) / (X[idx_dn_a]) +
                  2.0 * gdt * DN->Nsq() * dn_n_corr * Efdn * ABS(X[idx_dn_q]);

#if USEDIFF
  double dmdaup =
      lambda2 *
          (ff * OPT.Beta() * Qupsqdiva / (X[idx_up_a]) +
           OPT.Gravity() *
               (DN->GetDepth(X[idx_dn_a]) - UP->GetDepth(X[idx_up_a]) -
                (X[idx_up_a] + X[idx_dn_a]) * UP->GetDepthdA(X[idx_up_a]))) -
      gdt * UP->GetSR() + gdt * UP->Nsq() * up_n_corr * Qupsqabs * up_fric_da;

  double dmdadn =
      lambda2 *
          (-ff * OPT.Beta() * Qdnsqdiva / (X[idx_dn_a]) +
           OPT.Gravity() *
               (DN->GetDepth(X[idx_dn_a]) - UP->GetDepth(X[idx_up_a]) +
                (X[idx_up_a] + X[idx_dn_a]) * DN->GetDepthdA(X[idx_dn_a]))) -
      gdt * DN->GetSR() + gdt * DN->Nsq() * dn_n_corr * Qdnsqabs * dn_fric_da;
#else

  double dmdaup =
      lambda2 *
          (OPT.Beta() * Qupsqdiva / (X[idx_up_a]) +
           0.5 * OPT.Gravity() *
               (DN->GetDepth(X[idx_dn_a]) - UP->GetDepth(X[idx_up_a]) -
                (X[idx_up_a] + X[idx_dn_a]) * UP->GetDepthdA(X[idx_up_a]))) -
      gdt * UP->GetSR() + gdt * UP->Nsq() * up_n_corr * Qupsqabs * up_fric_da;

  double dmdadn =
      lambda2 *
          (-OPT.Beta() * Qdnsqdiva / (X[idx_dn_a]) +
           0.5 * OPT.Gravity() *
               (DN->GetDepth(X[idx_dn_a]) - UP->GetDepth(X[idx_up_a]) +
                (X[idx_up_a] + X[idx_dn_a]) * DN->GetDepthdA(X[idx_dn_a]))) -
      gdt * DN->GetSR() + gdt * DN->Nsq() * dn_n_corr * Qdnsqabs * dn_fric_da;
#endif

  M->CreateEntryNC(_row, UP->AIdx(), dmdaup);
  M->CreateEntryNC(_row, DN->AIdx(), dmdadn);
  M->CreateEntryNC(_row, UP->QIdx(), dmdqup);
  M->CreateEntryNC(_row, DN->QIdx(), dmdqdn);

  return OK;
}

// this function need to be synchronized with the evaluation code! it is
// becoming obsolete!!!
void MomStvEqn::PrintValues(FILE *F, double t, double dt, Node **n, double *X,
                            double *Xp, double *RHS) {
  double lambda2 = 2.0 * dt / _dx;
  double gdt = OPT.Gravity() * dt;
  Node *UP = n[_up_idx];
  Node *DN = n[_dn_idx];

  int idx_up_a = UP->AIdx();
  int idx_up_q = UP->QIdx();
  int idx_dn_a = DN->AIdx();
  int idx_dn_q = DN->QIdx();

  double Qupsq = X[idx_up_q] * X[idx_up_q];
  double Qdnsq = X[idx_dn_q] * X[idx_dn_q];
  double Qupsqdiva = Qupsq / X[idx_up_a];
  double Qdnsqdiva = Qdnsq / X[idx_dn_a];

  double Qupsqabs = X[idx_up_q] * ABS(X[idx_up_q]);
  double Qdnsqabs = X[idx_dn_q] * ABS(X[idx_dn_q]);

  double Efdn = DN->XS()->GetEqFriction(X[idx_dn_a]);
  double Efup = UP->XS()->GetEqFriction(X[idx_up_a]);

  /* !!!!!!
     need to worry about the depth correction term !!! */

  double dmdqup = 1.0 - 2.0 * lambda2 * OPT.Beta() * X[idx_up_q] / X[idx_up_a] +
                  2.0 * gdt * UP->Nsq() * Efup * ABS(X[idx_up_q]);
  double dmdqdn = 1.0 + 2.0 * lambda2 * OPT.Beta() * X[idx_dn_q] / X[idx_dn_a] +
                  2.0 * gdt * DN->Nsq() * Efdn * ABS(X[idx_dn_q]);

  double dmdaup =
      lambda2 * (OPT.Beta() * Qupsqdiva / X[idx_up_a] -
                 OPT.Gravity() * UP->XS()->GetCentroiddA(X[idx_up_a])) -
      gdt * UP->GetSR() +
      gdt * UP->Nsq() * Qupsqabs * UP->XS()->GetEqFrictiondA(X[idx_up_a]);

  double dmdadn =
      lambda2 * (-OPT.Beta() * Qdnsqdiva / X[idx_dn_a] +
                 OPT.Gravity() * DN->XS()->GetCentroiddA(X[idx_dn_a])) -
      gdt * DN->GetSR() +
      gdt * DN->Nsq() * Qdnsqabs * DN->XS()->GetEqFrictiondA(X[idx_dn_a]);

  fprintf(F, "** Qup=%6.2e Aup=%6.2e Qdown=%6.2e Adown=%6.2e\n", X[idx_up_q],
          X[idx_up_a], X[idx_dn_q], X[idx_dn_a]);
  fprintf(F,
          "MOMV: (%4d %4d)=%6.2e (%4d %4d)=%6.2e (%4d %4d)=%6.2e (%4d "
          "%4d)=%6.2e || (%4d)=%+6.2e\n",
          _row, UP->QIdx(), dmdqup, _row, DN->QIdx(), dmdqdn, _row, UP->AIdx(),
          dmdaup, _row, DN->AIdx(), dmdadn, _row, RHS[_row]);
}

void MomStvEqn::CreateEntry(SparseMatrix *w, Node **nn) {
  w->CreateEntry(_row, nn[_up_idx]->QIdx(), 0.0);
  w->CreateEntry(_row, nn[_up_idx]->AIdx(), 0.0);
  w->CreateEntry(_row, nn[_dn_idx]->QIdx(), 0.0);
  w->CreateEntry(_row, nn[_dn_idx]->AIdx(), 0.0);
}

// End
