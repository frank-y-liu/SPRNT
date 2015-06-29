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

/**********************************************************************
comments: the formulation should be taken to remove the known source 
    at the upstream boundary point to reduce the number of unknowns, which is 
    closely associated with the size of the matrix 
*******************************************************************/

#ifndef _SUBCATCH_H
#define _SUBCATCH_H

#include "node.h"
#include "sources.h"
#include "xsection.h"
#include "xstrans.h"
#include "ngraph.h"
#include "waveforms.h"
#include "solver_interface.h"

typedef int source_id;
typedef int node_id;

/**********************************************************************
  a "subcachment" is unit in which a system of SVE will be formulated
  and solved in one single matrix
 **********************************************************************/
class Subcatchment {
 protected:
  unsigned int     _cid;                 // id of the subcatchment
  int              _num_nodes;           // number of internal nodes
  int              _num_eqns;            // number of equations 
  int              _num_src;             // number of sources, can be Q, lat, or A

  int              _num_superc;          // number of super critical
  double           _max_froude;          // max of froude number
  double           _min_a;               // as the name says
  double           _min_q;               // 
  double           _max_a;               //
  double           _max_q;               //
  
  double           _t;                   // current time point
  double           _tm1;                 // previous time point
  double           _dt;                  // current time step
  
  GrowVec<Node*,     DRSS_NN_DFT>         _NS;    // collection of nodes;
  GrowVec<Equation*, DRSS_NN_DFT>         _EQ;    // collection of equations
  GrowVec<Source*,   DRSS_SRC_DFT>        _SRC;  // collection of Q source

  SparseMatrix                            _M;     // the jacobian matrix
  Solver_Interface                       *_s;

  /* internal storage, size varies */
  FlexVec<double>  _RHS;                 // the right-hand-side
  FlexVec<double>  _X;                   // solution			
  FlexVec<double>  _DX;                  // delta X, from 
  FlexVec<double>  _Xp;                  // solution @ previous newton step
  FlexVec<double>  _Xhigh;               // Upper limit
  FlexVec<double>  _Xlow;                // Lower limit;
  FlexVec<double>  _Xtm1;                // solution @ previous time point
  FlexVec<double>  _Froude;              // froude numbers
  FlexVec<double>  _Depth;               // to store the depth, computed from A
  FlexVec<double>  _Elevation;           // to store the elevation, computed from depth

  // temporary working space
  FlexVec<double,DRSS_NN_DFT>   _work;   // temporary working space
  FlexVec<double,DRSS_NN_DFT>   _workb;  // temporary working space

  int              _num_upsteam_bndy;    // number of upstream boundary conditions
  double           _down_a;              // downstream boundary conditions

  xs_trans         *_XST;                 // for translation
  ngraph            _G;                   // for topology check

  // private methods, hidden
  int NonLinearStep(double t, double dt, int jac_num, int max_iter, double tol, int do_bounding,
		    int &num_iter, double &minnorm, double &maxnorm, double &tdif, FILE *F=0);

  // evaluation methods
  // dt > 0 : unsteady
  // dt > 0 : steady, ignore entries related to Q.
  int EvaluateRHS(double t, double dt);
  int Evaluate(double t, double dt);

  /* utilities to compute a few others */
  int ComputeDepthAndElevation(int what2print);

 public:
  Subcatchment(int cid):_cid(cid),_num_nodes(0),_num_eqns(0),_num_src(0),
			_num_superc(0) { 
    _s=NULL;
    _max_froude=0.0; _t=0.0; _dt=0.0;
    _max_a = 0; _min_a = 0; _max_q=0; _min_q=0; _XST=NULL; }
  ~Subcatchment();

  /* get stats */
  int GetNumNodes() const { return _num_nodes;}
  int GetNumEqns() const { return _num_eqns; }
  int GetNumSrc() const { return _num_src; }
  int GetNumSuperC() const { return _num_superc; }
  double GetMaxFroude() const { return _max_froude; }

  void CalSolStat();
  double GetCurTime() const { return _t;  }
  double GEtCurDt() const { return _dt;   }
  double GetMaxA() const { return _max_a; }
  double GetMaxQ() const { return _max_q; }
  double GetMinA() const { return _min_a; }
  double GetMinQ() const { return _min_q; }

  /* expose some methods for graph */
  void InitGraph(int nn, int ns, int nj) { _G.Init(nn+5, ns+5, nj+5); }
  int  TopologyCheck() { return (_G.TChk() ); }
  void TopoPrintErrMsg(FILE *F, const char *h, char** table) {_G.PrintErrorMsg(F,h,table);}

  /* methods to construct the subcatchment, returns the node index*/
  node_id MakeNode(int id, float s0, float n, double x, double y, float q0, float a0, float z0, float h0, 
		   XsecType tp, double b0);  // RECT
  node_id MakeNode(int id, float s0, float n, double x, double y, float q0, float a0, float z0, float h0, 
		   XsecType tp, double b0, double s);  // TRAP
  node_id MakeNode(int id, float s0, float n, double x, double y, float q0, float a0, float z0, float h0, 
		   XsecType tp, int sec_num, double *xsecx, double *xsecy); // XY
  node_id MakeNode(int id, float s0, float n, double x, double y, float q0, float a0, float z0, float h0, 
		   XsecType tp, int num_pt, double *aa, double *pp, double *yy, double *ww); // INTRINSIC

  // determine the total number of nodes, and assign the indices to the X vector
  int AssignNodes();

  /* make time-varying sources, so that they can be later used to build equations
     since they are "dynamic", they can also be updated
     effectively only PWL sources are supported */
  // returns the index, return -1 if error
  source_id MakeSource(int sz, double *t, double *y, double t_scale=1.0, double y_scal=1.0); 
  // update, returns <0 if error
  int UpdateSource(source_id which, int sz, double *t, double *y, double t_scale=1.0, double y_scale=1.0);
  // retrieval 
  Source* GetSource(source_id which) { assert(which<_num_src); return _SRC[which];}


  // make equations, can be src, dep or stv
  int MakeStvEquation(int up_idx, int dn_idx, double dx, Source *s);
  int MakeArSrcEquation(int nidx, Source *s);
  int MakeQrSrcEquation(int nidx, Source *s);

  // overloaded version
  int MakeStvEquation(int up_idx, int dn_idx, double dx, source_id src=-1);
  int MakeArSrcEquation(int nidx, source_id src);
  int MakeQrSrcEquation(int nidx, source_id src);

  int MakeDepEquation(int nidx, int sz, int *indices, double *ratios);


  // Query methods
  Node* GetNode(int which) { assert(which<_num_nodes); return _NS[which]; }
  Equation* GetEquation(int which) { assert(which<_num_eqns); return _EQ[which]; }

  // method to allocate the solver
  int MakeSolver( int num_nodes );

  /* implementations need for Solve */
  double GetQNorm();   // norm of RHS, Q only
  double GetANorm();   // norm of RHS, A only
  double GetNorm();    // norm of RHS, both

  void GetDiffs(double &qdif, double &adif, double &tdif);  // diff between X and Xp
  int GetCelerity(double &max_cel, double &opt_dt);    // return celerity and best dt
  double CalFroude();                                  // return optimal step size and
						       // calculate froude numbers

  // as the name says
  void CopyXpToX(void) { memcpy(&_X[0], &_Xp[0], _num_eqns*sizeof(double) ); }
  void CopyXpToXtm1(void) { memcpy(&_Xtm1[0], &_Xp[0], _num_eqns*sizeof(double) ); }

  // bounding
  double Accept(double t, int needed, double& qdif, double &adif);
  int CheckMinA(double mina);   // check min A, not for bounding
  double QueryBreakPoint(double now, int &atbrkpt);
  
  /*  methods for solve */
  int InitSolutions(void);

  /* steady solve, this is the "psuedo unsteady" method */
  int SteadySolve(double dt, int jac_num, int max_iter, double tol );
  /* steady solve, when use_acc is set, use newton's method on A. 
     Usually it should be followed by a call to the psuedo unsteady method */
  int SteadySolve(int jac_num, int max_iter, double tol, int use_acc=1); 

  /* unsteady method
     if OUT != NULL : print results to OUT
     if WAVE != NULL : store the results to WAVE
   */
  int UnsteadySolve(double final_t, int jac_num, int max_iter, double tol, 
		    Waveforms *WV, FILE *OUT, int what2print, int pstart, char** N);

  /* printing */
  void PrintOnDemand(FILE *F, double tnow, int tstart, int what=0, char **NAMES=NULL); 

  /* methods for debugging */
  void PrintToFile(FILE *F, int converged_x=1);
  void PrintQToFile(FILE *F, int partial_list_only=1); // filter out those w/ zero xy
  void PrintAToFile(FILE *F, int partial_list_only=1);
  int SaveSteadyStateToFile(FILE *F);
  int LoadSteadyStateFromFile(FILE *F, int cmp_chksum=1);
  inline double GetQ(int nidx) { assert(nidx<_num_nodes); return (_Xp[2*nidx]);}
  inline double GetA(int nidx) { assert(nidx<_num_nodes); return (_Xp[2*nidx+1]); }
  void MatlabDumpDX(FILE *F, const char* name="DX");
  void MatlabDumpRHS(FILE *F, const char *name="RHS");
  void MatlabDumpX(FILE *F, const char *name="X");
  void MatlabDumpXp(FILE *F, const char *name="Xp");
  void PrintEvaluationValues(double t, double dt, FILE *F, const char *separator="NoName");
  void KinematicEstimate(FILE *F, char **N=NULL);
  int CheckBankFull();
  int CheckMinimalA();

  void PrintCoordToFile(FILE *F);
};

#endif
// Local Variables:
// mode: c++
// End:

