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

/********************************************************************
  definitions for St. Venant equations
  09/01/11      moved xsection to a separate file   
 ********************************************************************/

#ifndef _NODE_H
#define _NODE_H

#include <assert.h>

#include "xsection.h"
#include "flexvec.h"
#include "sources.h"
#include "matrix.h"

// Node definition
class Node {
 protected:
  /* non-critical info */
  unsigned int      _id;            // unique id
  double             _x;            // x-y coordiates in case we want to plot
  double             _y;            // consider change to int, need double precision
  double            _q0;            // initial flow rate
  double            _a0;            // inital wetted area

  /* we need these to work */
  float             _sr;            // bottom slope
  float             _nsq;           // manning's n, squared!!
  XSection         *_xs;            // cross section description, use pointers
  
  /* short cut to the index */
  int               _a_idx;         // index to the A unknown
  int               _q_idx;         // index to the Q unknown

  /* 
     definition of the sr line:
       1. The distance from datum to the bathymetry bottom is z0 (usually positive)
       2. The distance from datum to "slope" line is zR (usually positive)
       3. The distance from "slope" line to bathymetry bottom is hR (can be pos or neg)
             always true:  zR + hR = z0
       4. The distance from "slope" line to free surface is ha, this is computed by 
          solving SVE. Or y(A) + hR ( y(A) is computed from A, which is based on z0)
       5. The true depth is y = ha - hR, as reported
       6. The free surface elevation is zR + ha = z0 + y
       7. Sprnt will report y and free surface elevation (besides A)
       8. The slope of the "slope" line is sR, which should always be positive
   */
  float              _zr;           // used for reporting
  float              _hr;           // used to calculate depth

  // whether print
  int                _print;

  // debug flag
  static int        _debug;
  
 public:
 Node(unsigned int id, float sr, float n, XSection *xs):
  _id(id),_x(0.0),_y(0.0),_q0(0.0),_a0(0.0),_sr(sr),_xs(xs) {
    _nsq = n*n;
    _a_idx = -1;
    _q_idx = -1;
    _zr = 0.0;
    _hr = 0.0;
    _print = 0;
  };

  Node(unsigned int id, double sr, double n, XSection *xs,
       double x, double y, float q0, float a0, float zr, float hr, int p=1):
  _id(id),_x(x),_y(y),_q0(q0),_a0(a0),_sr(sr) {
    _xs = xs;
    _nsq = n*n;
    _q_idx = -1;
    _a_idx = -1;
    _zr = zr;
    _hr = hr;
    _print = p;
  }

  ~Node() { 
    if (_xs) delete _xs; 
    _xs = 0; 
  }

  /* access methods */
  int Id() const { return _id; }
  int QIdx() const { return _q_idx;}
  int AIdx() const { return _a_idx;}
  XSection* XS() const { return _xs;}
  double Nsq() const { return _nsq;}
  double GetSR() const { return _sr;}
  void SetSR(double sr) { assert(sr>0); _sr=sr;}
  double X() const { return _x;}
  double Y() const { return _y;}
  double& Q0() { return _q0;}
  double& A0() { return _a0;}
  float& HR()  { return _hr;}

  /* method to add more sutff */
  void AddLocationInfo(double x, double y);

  /* Evaluation methods */
  void AssignEntries(int q_idx, int a_idx);

  void SetInitValues(double *X);

  /* Kinematic routing based on the given value Q0 */
  double KinematicEstimate(double Q);

  /* print the froude number */
  double GetFroude(double *X);
  double GetFroude(double q, double a);

  // Depth is the solution of the SVE, which is relative to the sR line
  // AbsoluteDepth is the depth from the bathymetry bottom
  double GetDepth(double a) { return a < 0 ?-1.0 : _xs->GetDepth(a) + _hr; }
  double GetAbsoluteDepth(double a) { return a < 0 ? -1.0 : _xs->GetDepth(a); }
  double GetElevation(double a) { return _xs->GetDepth(a) + _zr + _hr; }
  double GetAbyDepth(double a) { return _xs->GetAbyDepth(a); }
  double GetDepthdA(double a) { return _xs->GetDepthdA(a); }

  // determine whether to print
  void SetPrint() { _print = 1; }
  void UnsetPrint() { _print = 0; }
  int Print() const { return _print;}

  static inline int Debug() { return _debug;}
  static inline void SetDebug() { _debug = 1; }
  static inline void ResetDebug() { _debug = 0; }
 
};

/* Equations */
class Equation {
 protected:
  int      _row;
  MMType   _type;

 public:
  Equation(int r):_row(r) {}
  virtual ~Equation() {}

  int Row() const { return _row;}
  MMType Type() const { return _type;}

  virtual void CreateEntry(SparseMatrix *w, Node** nodes) {}
  virtual int EvaluateRHS(double t, double dt, Node **n, double *X, double *Xp, double *RHS) { return 0; }
  virtual int EvaluateJac(double t, double dt, Node **n, double *X, double *Xp, SparseMatrix *M) { return 0; }
  virtual int Evaluate(double t, double dt, Node **n, double *X, double *Xp, double *RHS, SparseMatrix *M) { return 0; }

  virtual int EvaluateRHS(Node **n, double *X, double *Xp, double *RHS) { return 0; }
  virtual int Evaluate(Node **n, double *X, double *Xp, double *RHS, SparseMatrix *M) { return 0; }

  // for debugging only, could be quite slow
  virtual void PrintValues(FILE *F, double t, double dt, Node **n, double *X, double *Xp, double *RHS) {}
};

/**************************************
  independent A source, a.k.a. donwstream boundary conditions
 **************************************/
class ArSrcEqn : public Equation {
 private:
  int       _nidx;
  Source   *_S;
  static int _debug;

 public:
  ArSrcEqn(int r):Equation(r) {_type=A_SRC;}
  ~ArSrcEqn() { } // sources are owned by parent

  int EvaluateRHS(double t, double dt, Node **n, double *X, double *Xp, double *RHS);
  int EvaluateJac(double t, double dt, Node **n, double *X, double *Xp, SparseMatrix *M);
  int Evaluate(double t, double dt, Node **n, double *X, double *Xp, double *RHS, SparseMatrix *M);

  int EvaluateRHS(Node **n, double *X, double *Xp, double *RHS);
  int Evaluate(Node **n, double *X, double *Xp, double *RHS, SparseMatrix *M);

  void CreateEntry(SparseMatrix *w, Node** n);

  void Assign(int nidx) { assert(nidx>=0); _nidx = nidx; }
  void AddSource(Source *s) { _S = s;}

  void PrintValues(FILE *F, double t, double dt, Node **n, double *X, double *Xp, double *RHS);

  static inline int Debug() { return _debug;}
  static inline void SetDebug() { _debug = 1; }
  static inline void ResetDebug() { _debug = 0; }

};


/**************************************
  dependent A source, occurs at the junctions
 **************************************/
class ArDepEqn : public Equation {
 private:
  int                          _nidx;
  int                          _dep_idx;
  double                       _ratio;

 public:
  ArDepEqn(int r):Equation(r) {_type=A_DEP;}
  ~ArDepEqn() {};

  int EvaluateRHS(double t, double dt, Node **n, double *X, double *Xp, double *RHS);
  int EvaluateJac(double t, double dt, Node **n, double *X, double *Xp, SparseMatrix *M);
  int Evaluate(double t, double dt, Node **n, double *X, double *Xp, double *RHS, SparseMatrix *M);
  int EvaluateRHS(Node **n, double *X, double *Xp, double *RHS);
  int Evaluate(Node **n, double *X, double *Xp, double *RHS, SparseMatrix *M);

  void CreateEntry(SparseMatrix *w, Node** n);

  void Assign(int nidx) { assert(nidx>=0); _nidx = nidx; }  
  int AddDep(int dep_idx, double ratio);

  void PrintValues(FILE *F, double t, double dt, Node **n, double *X, double *Xp, double *RHS);
};

/**************************************
  mass part of SVE. a.k.a continuity equations
 **************************************/
class MasStvEqn : public Equation {
 private:
  int        _up_idx;
  int        _dn_idx;
  double     _dx;
  Source    *_pp;      // lateral flow, could be time varying, named pp
  static int _debug;
  static double _spin_coef;    // coefficient for spin-up

 public:
   MasStvEqn(int r):Equation(r) {_type=MAS_STV; _pp=NULL;}
  ~MasStvEqn() {  }   // sources are owned by parent

  int EvaluateRHS(double t, double dt, Node **n, double *X, double *Xp, double *RHS);
  int EvaluateJac(double t, double dt, Node **n, double *X, double *Xp, SparseMatrix *M);
  int Evaluate(double t, double dt, Node **n, double *X, double *Xp, double *RHS, SparseMatrix *M);
  int EvaluateRHS(Node **n, double *X, double *Xp, double *RHS);
  int Evaluate(Node **n, double *X, double *Xp, double *RHS, SparseMatrix *M);

  void CreateEntry(SparseMatrix *w, Node** n);

  void Assign(int upidx, int dnidx) { assert(dnidx>=0 && upidx>=0); _up_idx=upidx; _dn_idx=dnidx; }
  void AddDX(double dx) { _dx = dx;}
  double DX() const { return _dx;}
  int UpIdx() const { return _up_idx;}
  int DnIdx() const { return _dn_idx;}

  void AddLateralFlow(Source *s) { _pp = s; }
  Source *PP() { return _pp; }

  void PrintValues(FILE *F, double t, double dt, Node **n, double *X, double *Xp, double *RHS);

  static inline int Debug() { return _debug;}
  static inline void SetDebug() { _debug = 1; }
  static inline void ResetDebug() { _debug = 0; }
  
  static inline void SetSpinCoef(double c) { assert(c>=0.0 && c<=1.0); _spin_coef=c;}
  static inline double SpinCoef() { return _spin_coef; }

};

/**************************************
  independent Q source, a.k.a. upstream Q flow
 **************************************/
class QrSrcEqn : public Equation {
 private:
  int        _nidx;
  Source    *_S;
  static int _debug;

 public:
   QrSrcEqn(int r):Equation(r) {_type=Q_SRC;}
  ~QrSrcEqn() { } // source is owned by parent

  int EvaluateRHS(double t, double dt, Node **n, double *X, double *Xp, double *RHS);
  int EvaluateJac(double t, double dt, Node **n, double *X, double *Xp, SparseMatrix *M);
  int Evaluate(double t, double dt, Node **n, double *X, double *Xp, double *RHS, SparseMatrix *M);

  int EvaluateRHS(Node **n, double *X, double *Xp, double *RHS);
  int Evaluate(Node **n, double *X, double *Xp, double *RHS, SparseMatrix *M);


  void CreateEntry(SparseMatrix *w, Node** n);

  void Assign(int nidx) {assert(nidx>=0); _nidx = nidx;}
  void AddSource(Source *s) { _S = s;}
  Source* Src() { return _S; };

  void PrintValues(FILE *F, double t, double dt, Node **n, double *X, double *Xp, double *RHS);

  static inline int Debug() { return _debug;}
  static inline void SetDebug() { _debug = 1; }
  static inline void ResetDebug() { _debug = 0; }

};

/**************************************
  depdendent Q source, occurs at junctions
 **************************************/
class QrDepEqn : public Equation {
 private:
  int                          _nidx;
  int                          _num_deps;
  FlexVec<int, DRSS_DEP_DFT>   _deps;
  Source                      *_pp;       // pour port, could be NULL
  static int                   _debug;

 public:
   QrDepEqn(int r):Equation(r) { _type=Q_DEP; _pp=NULL; }
  ~QrDepEqn() {}   // source owned by parent

  int EvaluateRHS(double t, double dt, Node **n, double *X, double *Xp, double *RHS);
  int EvaluateJac(double t, double dt, Node **n, double *X, double *Xp, SparseMatrix *M);
  int Evaluate(double t, double dt, Node **n, double *X, double *Xp, double *RHS, SparseMatrix *M);

  int EvaluateRHS(Node **n, double *X, double *Xp, double *RHS);
  int Evaluate(Node **n, double *X, double *Xp, double *RHS, SparseMatrix *M);

  void CreateEntry(SparseMatrix *w, Node** n);

  void Assign(int nidx) { assert(nidx>=0); _nidx=nidx;}
  int AddDep(int sz, int* indices);
  
  void AddPourPort(Source *s) { _pp = s;}
  Source* PP() { return _pp;}

  void PrintValues(FILE *F, double t, double dt, Node **n, double *X, double *Xp, double *RHS);
  static inline int Debug() { return _debug;}
  static inline void SetDebug() { _debug = 1; }
  static inline void ResetDebug() { _debug = 0; }

};

/**************************************
  momentum portion of SVE. a.k.a. dynamic equation
 **************************************/
class MomStvEqn : public Equation {
 private:
  int          _up_idx;
  int          _dn_idx;
  double       _dx;

 public:
   MomStvEqn(int r):Equation(r) { _type = MOM_STV;}
  ~MomStvEqn() {}

  int EvaluateRHS(double t, double dt, Node **n, double *X, double *Xp, double *RHS);
  int EvaluateJac(double t, double dt, Node **n, double *X, double *Xp, SparseMatrix *M);
  int Evaluate(double t, double dt, Node **n, double *X, double *Xp, double *RHS, SparseMatrix *M);

  int Evaluate(Node **n, double *X, double *Xp, double *RHS, SparseMatrix *M);
  int EvaluateRHS(Node **n, double *X, double *Xp, double *RHS);

  void CreateEntry(SparseMatrix *w, Node** n);

  void Assign(int up, int dn) { assert(up>=0 && dn>=0); _up_idx=up; _dn_idx=dn;}
  void AddDX(double dx) { _dx = dx;}  

  void PrintValues(FILE *F, double t, double dt, Node **n, double *X, double *Xp, double *RHS);
};

#endif

// Local Variables:
// mode: c++
// End:

