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

/******************************************************************************
  methods related to the subcathment class. here is where the system is explicitly
  solved
******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "subcatch.h"
#include "node.h"
#include "xsmodel.h"
#include "blaspr.h"

#ifndef ISNAN
#define ISNAN( A )  ( !( (A)==(A) ) )
#endif

#define PRT_FMT " %.6e"
#define PRT_XY_FMT " %.4f %.4f"

Subcatchment::~Subcatchment() {
  if (_num_nodes > 0 ) {
    for (int j=0; j<_num_nodes; j++) {
      if (_NS[j] > 0 ) delete _NS[j];
      _NS[j] = 0;
    }
  }

  if (_num_eqns > 0 ) {
    for (int j=0; j<_num_eqns; j++) {
      if ( _EQ[j] > 0 ) delete _EQ[j];
      _EQ[j] = 0;
    }
  }

  if ( _num_src > 0 ) {
    for (int jj=0; jj<_num_src; jj++) {
      if ( _SRC[jj] ) delete _SRC[jj];
      _SRC[jj]=0;
    }
  }

  _num_nodes = 0;
  _num_eqns = 0;

  if ( _s ) delete _s;
  _s = 0;

  if ( _XST ) delete _XST;
  _XST = 0;
}

/* as the name says, print them  as:

   index x-coordinate y-coordinate

*/
void Subcatchment::PrintCoordToFile(FILE *F) {
  int jj;

  for (jj=0; jj<_num_nodes; jj++) {
    fprintf(F, "%d %.4f %.4f\n", jj, _NS[jj]->X(), _NS[jj]->Y() );
  }
}

void Subcatchment::PrintToFile(FILE *F, int previous_x) {
  fprintf(F, "%%%%%% index  Q   A\n");

  if (previous_x == 1) {
    for (int i=0; i<_num_eqns; i+=2) fprintf(F, "%5d  %10.5e %10.5e\n", i/2,_X[i],_X[i+1]);
  } else {
    for (int i=0; i<_num_eqns; i+=2) fprintf(F, "%5d  %10.5e %10.5e\n", i/2,_Xp[i],_Xp[i+1]);
  }
}

/* if partial_list_only==1, only print those with nonzero xy coordiates */
void Subcatchment::PrintQToFile(FILE *F, int partial_list_only) {
  int j;
  if (partial_list_only==1) {
    for (j=0; j<_num_nodes; j++) {
      if ( _NS[j]->X() > 0 || _NS[j]->Y() > 0 ) fprintf(F, "%.5e\n", _X[2*j]);
    }
  } else {
    for (j=0; j<_num_nodes; j++) fprintf(F, "%.5e\n", _Xp[2*j]);
  }
}

void Subcatchment::PrintAToFile(FILE *F, int partial_list_only) {
  int j;
  if (partial_list_only==1) {
    for (j=0; j<_num_nodes; j++) {
      if ( _NS[j]->X() > 0 || _NS[j]->Y() > 0 ) fprintf(F, "%.5e\n", _X[2*j+1]);
    }
  } else {
    for (j=0; j<_num_nodes; j++) fprintf(F, "%.5e\n", _Xp[2*j+1]);
  }
}

/* save steady state to a file 
   save Xtm1 - for steady state, they are the same as X;
   but it also has the benefit of capable of having a hot restart for unsteady case
 */
int Subcatchment::SaveSteadyStateToFile(FILE *F) {
  int j;
  fprintf(F,"%d %s\n", _num_eqns, STAT.ChkSum() );
  for (j=0; j<_num_eqns; j++) fprintf(F,"%.29e\n", _Xtm1[j]);
  return OK;
}

/* load the steady state from file, overwrite the existing values
   a quick check to see if the sizes match */
int Subcatchment::LoadSteadyStateFromFile(FILE *F, int cmp_chksum) {
  int sz, j;
  char buf[512], chksum[10];
  double v;

  if ( (fgets(buf,512,F)) == NULL ) {
    fprintf(stderr,"[WW]: Unable to read content from steadystate file %s. Stored steady state ignored!\n", STAT.SSFile() );
    return WARNING;
  }

  sscanf(buf, "%d %s", &sz, chksum);

  if ( sz != _num_eqns ) {
    fprintf(stderr,"[WW]: Size mismatch in steadystate file %s. Stored steady state ignored!\n", STAT.SSFile() );
    return WARNING;
  }

  if ( cmp_chksum==1) {
    if ( strncmp(chksum, STAT.ChkSum(), 8) != 0 ) {
      fprintf(stderr,"[WW]: Signature mismatch found in steadystate file %s. Stored steady state ignored!\n", STAT.SSFile() );
      return WARNING;
    }
  } // else we skip it

  j=0;
  while ( (fgets(buf,512,F)) != NULL ) {
    v = atof(buf);
    if ( isnan(v) || isinf(v) || ( j%2==1 && v<1e-8 ) ) {
      fprintf(stderr,"[WW]: Found bad values in ssfiles. Ignored\n");
      return WARNING;
    }
    _Xp[j++] = v;
  }

  if ( j!=_num_eqns) {
    fprintf(stderr,"[WW]: Not sufficient values in ssfile. Ignored\n");
    return WARNING;
  }
  return OK;
}

void Subcatchment::MatlabDumpDX(FILE *F, const char *name) {
  int j;
  fprintf(F,"%s = [\n", name);
  for (j=0; j<2*_num_nodes; j++) fprintf(F,"%.4e\n", _DX[j]);
  fprintf(F,"];\n");
}

void Subcatchment::MatlabDumpX(FILE *F, const char *name) {
  int j;
  fprintf(F,"%s = [\n", name);
  for (j=0; j<2*_num_nodes; j++) fprintf(F,"%.4e\n", _X[j]);
  fprintf(F,"];\n");
}

void Subcatchment::MatlabDumpXp(FILE *F, const char *name) {
  int j;
  fprintf(F,"%s = [\n", name);
  for (j=0; j<2*_num_nodes; j++) fprintf(F,"%.4e\n", _Xp[j]);
  fprintf(F,"];\n");
}

void Subcatchment::MatlabDumpRHS(FILE *F, const char *name) {
  int j;
  fprintf(F,"%s = [\n", name);
  for (j=0; j<2*_num_nodes; j++) fprintf(F,"%.4e\n", _RHS[j]);
  fprintf(F,"];\n");
}

void Subcatchment::CalSolStat() {
  int j;
  _min_a = 9999.99;
  _max_a = -9999.99;
  _min_q = 99999.99;
  _max_q = -9999.99;

  for (j=0; j<_num_eqns; j+=2 ) {
    _min_q = _min_q < _X[j] ? _min_q : _X[j];
    _max_q = _max_q > _X[j] ? _max_q : _X[j];
  }

  for (j=1; j<_num_eqns; j+=2 ) {
    _min_a = _min_a < _X[j] ? _min_a : _X[j];
    _max_a = _max_a > _X[j] ? _max_a : _X[j];
  }

}
// Make a node with the given xsection description, for RECT type only
//  when use for RECT, b0 is the bottom width, 
// returns node index if successful, -1 if failed ( we will never fail )
node_id Subcatchment::MakeNode(int id, float s0, float n, double x, double y, float q0, float a0, float z0, float h0, 
			   XsecType tp, double b0) {
  R_XSection *rx;
  Node  *node;

  assert ( tp == RECT);
  assert ( b0 > 0 );

  rx = new R_XSection(0, b0);
  node = new Node(id, s0, n,(XSection*)rx, x, y, q0, a0, z0,h0);

  _NS[_num_nodes++] = node;

  _G.AddNode( id ); // insert to graph
  return (_num_nodes-1);
}

// Make a node with the given xsection description, for TRAP type only
node_id Subcatchment::MakeNode(int id, float s0, float n, double x, double y, float q0, float a0, float z0, float h0,XsecType tp, double b0, double s) {
  T_XSection *tx;
  Node  *node;

  assert ( tp == TRAP );
  assert ( b0 > 0 );

  tx = new T_XSection(0, b0,s);
  node = new Node(id, s0, n, (XSection*)tx, x, y, q0, a0, z0, h0);

  _NS[_num_nodes++] = node;

  _G.AddNode( id ); // insert to graph

  return (_num_nodes-1);
}

// make a node for the XY xsection
node_id Subcatchment::MakeNode(int id, float s0, float n, double x, double y, float q0, float a0, float z0, float h0,
			   XsecType tp, int sec_num, double *xsecx, double *xsecy) {
  XY_XSection *tx;
  Node        *node;
  char        buf[512];
  assert ( tp == SPLINE );
  assert ( sec_num > 2);

  sprintf(buf,"node %d :", id);
  OPT.CopyToBuf(buf);

  if ( !_XST ) _XST = new xs_trans;

  tx = new XY_XSection(id);
  tx->Build(sec_num, xsecx, xsecy, _XST);
  node = new Node(id, s0, n, (XSection*)tx, x,y,q0, a0, z0, h0);

  // printf("node %d, bf=%.3f\n", id, node->XS()->GetBkfA());
  _NS[_num_nodes++] = node;
  _G.AddNode( id ); // insert to graph
  return (_num_nodes-1);
}

// make a node for the INTRINSIC xsection
// the underline type is actuall XY!x
node_id Subcatchment::MakeNode(int id, float s0, float n, double x, double y, float q0, float a0, float z0, float h0,
			       XsecType tp, int num_pts, double *aa, double *pp, double *yy, double *ww) {

  XY_XSection *tx;
  Node        *node;
  char        buf[512];
  assert ( tp == INTRINSIC );
  assert ( num_pts > 5 );

  sprintf(buf,"node %d :", id);
  OPT.CopyToBuf(buf);

  tx = new XY_XSection(id);
  // build the cross section
  tx->Build(num_pts, aa, pp, yy, ww);
  node = new Node(id, s0, n, (XSection*)tx, x,y, q0, a0, z0, h0);

  // printf("node %d, bf=%.3f\n", id, node->XS()->GetBkfA());
  _NS[_num_nodes++] = node;
  _G.AddNode( id ); // insert to graph
  return (_num_nodes-1);
}

// assign the node index to the nodes
int Subcatchment::AssignNodes() {
  for (int j=0; j<_num_nodes; j++) {
    _NS[j]->AssignEntries(2*j, 2*j+1);  // the 1st is Q, the 2nd is A
  }
  return OK;
}

/* assemble the time-varying sources so that they can be later used 
   either as Q term, lateral term, or boundary A
   So far only PWL is supported
 */
source_id Subcatchment::MakeSource(int sz, double *t, double *y, double t_scale, double y_scale) {
  PWL *ns = new PWL;
  int rc;
  for (int jj=0; jj<sz; jj++) {
    t[jj] *= t_scale;
    y[jj] *= y_scale;
  }
  rc = ns->Assign(sz, t, y);
  if ( rc < 0 ) return (-1); 

  _SRC[_num_src] = (Source*) ns;
  return (_num_src++);
}

// rc < 0 if error
int Subcatchment::UpdateSource(source_id which, int sz, double *t, double *y, double t_scale, double y_scale){
  assert ( which < _num_src );

  if ( _t > t[0] ) return (-1);  // we only allow update for the future
  if ( PWL::CheckTime(sz, t) < 0 ) return (-1);
  
  for (int jj=0; jj<sz; jj++) {
    t[jj] *= t_scale;
    y[jj] *= y_scale;
  }

  PWL *ns = (PWL*) _SRC[which];

  return (ns->Assign(sz, t, y));
}

// make equations for the independent sources
// returns the number of equations added
int Subcatchment::MakeArSrcEquation(int nidx, Source *s) {
  ArSrcEqn *E = new ArSrcEqn(_num_eqns);
  E->Assign(nidx);
  E->AddSource(s);
  
  _EQ[_num_eqns++] = E;
  
  s->SetParent( nidx ); // so nobody can used it again
  _G.InsertA(nidx, s->Evaluate(0));
  return (1);
}

int Subcatchment::MakeArSrcEquation(int nidx, source_id src) {
  Source* p=this->GetSource(src);
  return (MakeArSrcEquation(nidx, p));
}

int Subcatchment::MakeQrSrcEquation(int nidx, Source *s) {
  QrSrcEqn *E = new QrSrcEqn(_num_eqns);
  E->Assign(nidx);
  E->AddSource(s);
  _EQ[_num_eqns++] = E;
  s->SetParent( nidx ); // prevent from mis-use
  _G.InsertQ(nidx, s->Evaluate(0) );
  return (1);
}

int Subcatchment::MakeQrSrcEquation(int nidx, source_id src) {
  Source *p = this->GetSource(src);
  return (MakeQrSrcEquation(nidx, p));
}

// make equation for the dependent sources
// make sure we have do it for both upstream and downstream nodes
// returns the number of equations added
int Subcatchment::MakeDepEquation(int dn_idx, int up_sz, int *up_idx, double *ratios) {
  // first from up stream to down stream
  for (int j=0; j<up_sz; j++) {
    ArDepEqn *E = new ArDepEqn(_num_eqns);
    E->Assign(up_idx[j]);
    E->AddDep(dn_idx, ratios[j]);
    _EQ[_num_eqns++] = E;
  }
  
  QrDepEqn *E = new QrDepEqn(_num_eqns);
  E->Assign(dn_idx); 
  E->AddDep(up_sz, up_idx);
  _EQ[_num_eqns++] = E;

  _G.AddJunction(up_sz, up_idx, ratios, dn_idx);
  return (up_sz+1);
}

// add st venant equations, need to call it twice, one for momentum and one for mass
// returns the number of equation added (e.g., 1 means good )
int Subcatchment::MakeStvEquation(int up_idx, int dn_idx, double dx, Source *pp) {
  MasStvEqn *E1;
  MomStvEqn *E2;

  assert ( dx > 0.0 );

  E1 = new MasStvEqn(_num_eqns);
  E1->Assign(up_idx, dn_idx);
  E1->AddDX(dx);
  _EQ[_num_eqns++] = E1;
  
  E2 = new MomStvEqn(_num_eqns);
  E2->Assign(up_idx, dn_idx);
  E2->AddDX(dx);
  _EQ[_num_eqns++] = E2;
  
  if ( pp ) { 
    E1->AddLateralFlow( pp );  // only the mass eqn is affected by pp
    pp->SetParent( up_idx );   // this would prevent the src to be used by anybody else

    // we only insert lateral flow if we do NOT use spin up
    if ( OPT.SpinUpTime() <= 0 )  _G.InsertQ(up_idx, pp->Evaluate(0));
  }

  _G.AddEdge(up_idx, dn_idx, dx);
  return (2);
}

int Subcatchment::MakeStvEquation(int up_idx, int dn_idx, double dx, source_id src) {
  Source *p = (src==-1)? NULL : this->GetSource(src);
  return ( MakeStvEquation(up_idx, dn_idx, dx, p));
}

// allocate the solver, and all others
int Subcatchment::MakeSolver(int mx_sz) {
  int j;
  Node **nn = (Node**)_NS;

  if ( _num_eqns - 2*_num_nodes != 0 ) {
    printf("The numbers of equations and unknowns are not the same (#eqn=%d, #unk=%d). Bummer!\n", 
	   _num_eqns, 2*_num_nodes);
    return ERROR;
  }
  
  // create entries in the matrix
  _M.MakeSpace(5*_num_eqns);    // sparsity is 4 per entry, give or take a few
  for (j=0; j<_num_eqns; j++) {
    _EQ[j]->CreateEntry(&_M, nn);
  }

  // allocate the vectors
  _Froude.size(_num_nodes);
  _Depth.size(_num_nodes);
  _Elevation.size(_num_nodes);

  _RHS.size(_num_eqns);
  _X.size(_num_eqns);
  _DX.size(_num_eqns);
  _Xp.size(_num_eqns);
  _Xhigh.size(_num_eqns);
  _Xlow.size(_num_eqns);
  _Xtm1.size(_num_eqns);

  _M.SetDimensions(_num_eqns, _num_eqns);

  _s = new Solver_Interface;
  // heuristics to allocate different solvers
#if 1 
  if ( mx_sz < 16383 ) {
    if ((_s->Setup(L_KLU)) < 0 ) {
      printf("cannot make solver!\n");
      return (-1);
    }
  } else {
    if ((_s->Setup(L_UMF)) < 0 ) return (-1);
  }
#else
  if ((_s->Setup(L_KLU))< 0 ) return (-1);
#endif

  _s->Init(_M.Nnz(), _M.Ndim(0), _M.Row(), _M.Col(), _M.Data());
  printf("[II]: successfully allocated solver, type = %d, size = %d\n", _s->SolverType(),_num_eqns );

  return OK;
}

double Subcatchment::GetQNorm() {
  // we implicitly assume Q start from 0, A start from 1
  // this code has to be updated if this assumption is broken
  const int twoi = 2;
  int n = _num_eqns;
  
  return ( FORTRAN(dnrm2)(&n, &_RHS[0], &twoi ) );
}

double Subcatchment::GetANorm() {
  // we implicitly assume Q start from 0, A start from 1
  // this code has to be updated if this assumption is broken
  const int twoi = 2;
  int n = _num_eqns;
  
  return ( FORTRAN(dnrm2)(&n, &_RHS[1], &twoi ) );

}

double Subcatchment::GetNorm() {
  const int onei = 1;
  int n = _num_eqns;
  
  return ( FORTRAN(dnrm2)(&n, &_RHS[0], &onei ) );
}


void Subcatchment::GetDiffs(double &qdiff, double &adiff, double &tdiff) {
  const double negone = -1.0;
  const int onei = 1;
  const int twoi = 2;

  _work.size(_num_eqns);
  memcpy(&_work[0], &_X[0], _num_eqns*sizeof(double) );
  FORTRAN(daxpy)(&_num_eqns, &negone, &_Xtm1[0], &onei, &_work[0], &onei);
  qdiff = FORTRAN(dnrm2)(&_num_eqns, &_work[0], &twoi);
  adiff = FORTRAN(dnrm2)(&_num_eqns, &_work[1], &twoi);
  tdiff = sqrt( qdiff*qdiff + adiff*adiff);
}

// calculate froude number, also returns the ideal step size based on the average of 
// flow velocity and celerity
double Subcatchment::CalFroude() {
  int j, cntr, bl;
  const int onei = 1;
  const double sqg = OPT.SqrtG();
  MasStvEqn *ME;
  double tiny, v;

  _work.size(_num_eqns);
  _workb.size(_num_nodes);
  
  // _work holds the velociyt at even locations
  memcpy(&_work[0], &_X[0], _num_eqns*sizeof(double) );
  for (j=0; j<_num_eqns; j+=2 ) _work[j] /= _work[j+1]; 
  
  // _workb holds celerity, sqrt(g*D), where D is the hydraulic depth
  for (j=0; j<_num_nodes; j++ ) _workb[j] = _X[_NS[j]->AIdx()] / _NS[j]->XS()->GetWidth( _X[_NS[j]->AIdx()] );
  for (j=0; j<_num_nodes; j++ ) _workb[j] = sqrt(_workb[j]);
  FORTRAN(dscal)(&_num_nodes, &sqg, &_workb[0], &onei);
  
  // calculate the froude number
  for (j=0; j<_num_nodes; j++) _Froude[j] = _work[2*j]/_workb[j];

  // find out the max froude number, and how many of them are super critical, 
  // could be expensive
  tiny = -1.0;
  cntr = 0;
  bl = -1;
  for (j=0; j<_num_nodes; j++) {
    if ( _Froude[j] > tiny        ) { bl = j; tiny = _Froude[j]; }
    if ( _Froude[j] > OPT.SuperC()) cntr++;
  }
  _num_superc = cntr;
  _max_froude = _Froude[bl];

  // calculate step size, use fudge number in OPT to scale it up
  // we use a heuristic method in place of the Froude number
  int tmp;
  tiny = 1e49;
  for (j=0; j<_num_eqns; j++) {
    if ( _EQ[j]->Type() != MAS_STV ) continue;
    ME = (MasStvEqn*)_EQ[j];
    tmp = ME->UpIdx();
    v = ME->DX()/(_workb[tmp] + _work[2*tmp]);
    tiny = (v < tiny && v>0) ? v : tiny;
  }

  return ( OPT.MaxDtR()*tiny );   // magic number

}

// calculate celerity and optimal dt, we actually need the froude number!
int Subcatchment::GetCelerity(double &max_cel, double &op_dt) {
  int cntr;
  int j, bl;
  const int onei = 1;
  const double sqg = OPT.SqrtG();
  MasStvEqn *ME;
  double tiny, v;
  int n=_num_nodes-1;

  _work.size(_num_eqns);
  _workb.size(_num_nodes);

  // _work holds the velociyt at even locations
  memcpy(&_work[0], &_X[0], _num_eqns*sizeof(double) );
  for (j=0; j<_num_eqns; j+=2 ) _work[j] /= _work[j+1]; 

  // _workb holds celerity, sqrt(g*D), where D is the hydraulic depth
  for (j=0; j<_num_nodes; j++ ) _workb[j] = _X[_NS[j]->AIdx()] / _NS[j]->XS()->GetWidth( _X[_NS[j]->AIdx()] );
  for (j=0; j<_num_nodes; j++ ) _workb[j] = sqrt(_workb[j]);
  FORTRAN(dscal)(&_num_nodes, &sqg, &_workb[0], &onei);
  
  // calculate the froude number
  for (j=0; j<_num_nodes; j++) _Froude[j] = _work[2*j]/_workb[j];

  // find out the max froude number, and how many of them are super critical, 
  // could be expensive
  tiny = -1.0;
  cntr = 0;
  for (j=0; j<_num_nodes; j++) {
    if ( _Froude[j] > tiny ) { bl = j; tiny = _Froude[j]; }
    if (_Froude[j] > OPT.SuperC()) cntr++;
  }
  if ( j < _num_nodes ) {
    _num_superc = cntr;
    _max_froude = _Froude[bl];
  }

  // find the max celerity, reuse bl
  bl = FORTRAN(idamax)(&n, &_workb[0], &onei);

  int tmp;
  tiny = 1e49;
  for (j=0; j<_num_eqns; j++) {
    if ( _EQ[j]->Type() != MAS_STV ) continue;
    ME = (MasStvEqn*)_EQ[j];
    tmp = ME->UpIdx();
    v = ME->DX()/(_workb[tmp] + _work[2*tmp]);
    tiny = (v < tiny && v>0) ? v : tiny;
  }

  max_cel = _workb[bl];
  op_dt = tiny;

  return (max_cel>OPT.MaxC() ? 1 : 0);
}

// Kinematic routing based on Q0, make sure they are populated before calling this
// function
// whatever in A will be ignored
void Subcatchment::KinematicEstimate(FILE *F, char** NAMES) {
  int jj;
  double a,y,fn,q;
  if ( OPT.UseMetric() == 1 ) {
    fprintf(F,"%%%% %6s %9s %9s %9s %9s\n", "id", "Q0(m3/s)", "A0(m2)", "Depth(m)", "Froude Num");
    for (jj=0; jj<_num_nodes; jj++) {
      q = _NS[jj]->Q0();
      a = _NS[jj]->KinematicEstimate( q );
      y = _NS[jj]->GetDepth(a);
      fn = _NS[jj]->GetFroude(q, a);
      if (!NAMES) fprintf(F,"%9d ", _NS[jj]->Id());
      else        fprintf(F,"%9s ", NAMES[ _NS[jj]->Id() ] );
      fprintf(F,"%9.3f", _NS[jj]->Q0());
      if ( a > 0 ) {
	fprintf(F," %9.3f %9.3f %9.3f\n", a, y, fn);
      } else {
	fprintf(F," %9s %9s %9s\n", "N/A", "N/A", "N/A");
      }
    }
  } else {
    fprintf(F,"%%%% %6s %9s %9s %9s %9s\n", "id", "Q0(f3/s)", "A0(f2)", "Depth(f)", "Froude Num");
    for (jj=0; jj<_num_nodes; jj++) {
      q = _NS[jj]->Q0();
      a = _NS[jj]->KinematicEstimate( q );
      y = _NS[jj]->GetDepth(a);
      fn = _NS[jj]->GetFroude(q, a);
      if (!NAMES) fprintf(F,"%9d ", _NS[jj]->Id());
      else        fprintf(F,"%9s ", NAMES[ _NS[jj]->Id() ] );
      fprintf(F,"%9.3f", _NS[jj]->Q0()/OPT.F3toM3());
      if ( a > 0 ) {
	fprintf(F," %9.3f %9.3f %9.3f\n", a/OPT.F2toM2(), y/OPT.FtoM(), fn);
      } else {
	fprintf(F," %9s %9s %9s\n", "N/A", "N/A", "N/A");
      }
    }
  }
}

// initialize the values
int Subcatchment::InitSolutions(void) {
  // each node knows how to set itself up, note we set to _Xp
  for (int j=0; j<_num_nodes; j++) _NS[j]->SetInitValues((double*)_Xp);

  // also copy to Xtm1 and reset RHS
  memcpy((double*)_Xtm1, (double*)_Xp, sizeof(double)*2*_num_nodes);
  memset((double*)_RHS, 0, sizeof(double)*2*_num_nodes);
  
  return OK;
}

void Subcatchment::PrintEvaluationValues(double t, double dt, FILE *F, const char *name) {
  int j;
  for (j=0; j<_num_eqns; j++) {
    _EQ[j]->PrintValues(F, t, dt, (Node**)_NS, (double*)_X, (double*)_Xp, (double*)_RHS);
  }
}

// we use dt to distinguish whether evaluate is for steady and unsteady
//  for unsteady solve - both Q and A entries are valid
//  for steady solve - only A entries are valid, the Q entries are not solved
//                     they are taken from graph traversal
int Subcatchment::Evaluate(double t, double dt) {
  _M.Clear();

  _s->Clear();

  if ( dt > 0 ) {  
    for (int j=0; j<_num_eqns; j++) {
      _EQ[j]->Evaluate(t, dt, (Node**)_NS, (double*)_X, (double*)_Xp, (double*)_RHS, &_M);
    }
  } else {
    for (int j=0; j<_num_eqns; j++) {
      _EQ[j]->Evaluate((Node**)_NS, (double*)_X, (double*)_Xp, (double*)_RHS, &_M);
    }
  }
  _M.SetDimensions(_num_eqns, _num_eqns);
  _s->Init(_M.Nnz(), _M.Ndim(0), _M.Row(), _M.Col(), _M.Data());
  return OK;
}

int Subcatchment::EvaluateRHS(double t, double dt) {
  if ( dt > 0 ) {
    for (int j=0; j<_num_eqns; j++) {
      _EQ[j]->EvaluateRHS(t, dt, (Node**)_NS, (double*)_X, (double*)_Xp, (double*)_RHS);
    }
  } else {
    for (int j=0; j<_num_eqns; j++) {
      _EQ[j]->EvaluateRHS((Node**)_NS, (double*)_X, (double*)_Xp, (double*)_RHS);
    }
  }
  _M.SetDimensions(_num_eqns, _num_eqns);
  return OK;
}

// calculate the qdiff, adiff and tdiff, returns tdiff
// does the copying and updates the time,
// also returns difference between the previous solutions if needed
// if needed == 0, accept only
// if needed == 1, returns the relative difference
// if needed == 2, returns total difference + qdifference + adifferene
double Subcatchment::Accept(double t, int needed, double& qdiff, double &adiff) {
  const int twoi=2;
  const double negone = -1.0;
  const int onei = 1;
  double tdiff=OPT.BogusLimit();
  
  qdiff = OPT.BogusLimit();
  adiff = OPT.BogusLimit();

  if ( needed > 0 ) {
    // relative norm
    _work.size(_num_eqns);
    memcpy(&_work[0], &_X[0], _num_eqns*sizeof(double) );
    FORTRAN(daxpy)(&_num_eqns, &negone, &_Xp[0], &onei, &_work[0], &onei);

    // relative difference, only meaningful when all positive 
    for (int j=0; j<_num_eqns; j++) _work[j] /= _Xp[j];  

    if (needed >= 2) {
      qdiff = FORTRAN(dnrm2)(&_num_nodes, &_work[0], &twoi);  // skip 2, the number should be half!
      adiff = FORTRAN(dnrm2)(&_num_nodes, &_work[1], &twoi);
      tdiff = sqrt( qdiff*qdiff + adiff*adiff);
    } else {
      tdiff = FORTRAN(dnrm2)(&_num_eqns, &_work[0], &onei);
    }
  }

  // copy and accept the time point
  memcpy(&_Xtm1[0], &_Xp[0], _num_eqns*sizeof(double) ); // save a copy for printing
  _tm1 = _t;
  memcpy(&_Xp[0], &_X[0], _num_eqns*sizeof(double) );
  _t = t;

  // print the results for debugging, 
  // WARNING: there could be a LOT of printout. Better pipe to a file
  if ( OPT.DebugLevel() >= 2 ) {
    double qq, aa, eqf, sf;
    fprintf(stderr,"   T = %.5e\n", _t);
    for (int j=0; j<_num_nodes; j++) {
      qq = _X[_NS[j]->QIdx()];
      aa = _X[_NS[j]->AIdx()];
      eqf = _NS[j]->XS()->GetEqFriction(aa);
      sf = qq*qq*_NS[j]->Nsq()*eqf/aa;
      fprintf(stderr,"      %5d Q( %11.4e ) A( %11.4e ) Y( %11.4e ) EqFrc( %11.4e ) s0( %11.4e ) sf( %11.4e )\n",
	     _NS[j]->Id(), qq, aa,  _NS[j]->GetDepth( aa ),
	     eqf, _NS[j]->GetSR(), sf);
    }
  }
  return tdiff;
}

// check whether A is less than mina
int Subcatchment::CheckMinA(double mina) {
  int rc = OK;

  for (int j=1; j<_num_eqns; j+=2) {
    if (_X[j]-mina < 0 ) {
#if 0
      _X[j] = -0.1*_X[j];
#else
      rc = ERROR; 
      break;
#endif
    }
  }
  return rc;
}

// simple linear interpolation, nothing sophisticated
// returns evaluation of xx, assuming we have (t1,x1) and (t2,x2) 
inline double lin_interp(double t1, double x1, double t2, double x2, double tt) {
  return ( x1 + (tt-t1)*(x2-x1)/(t2-t1) );
}

// what to print:
//    four bits 
//    0x08         Q
//    0x04         A
//    0x02         D
//    0x01         Z
// combination are allowed, but the sequence is fixed
void Subcatchment::PrintOnDemand(sptFile F, double tnow, int pstart, int what, char **NAMES) {
  double fq, fa, fd;
  if (OPT.UseMetric() ==1 ) {
    fq = 1.0;
    fa = 1.0;
    fd = 1.0;
  } else {
    fq = 1. / OPT.F3toM3();
    fa = 1. / OPT.F2toM2();
    fd = 1. / OPT.FtoM();
  }

  static int prev_t = pstart-OPT.PrintInterval();

  if ( ! (what & ( PRT_Q | PRT_A | PRT_D | PRT_Z)) ) return;  

  int interval = OPT.PrintInterval();
  int tnowmin = (int)(tnow/60.0);

  if (interval == 0 ) {  // we print whatever we have if interval is 0
    // print the header
    if (OPT.UseMetric() ==1 ) {
      sptFprintf(F,"*** id time(min)");
      if ( what & PRT_Q ) sptFprintf(F," flow(m3/s)");
      if ( what & PRT_A ) sptFprintf(F," wet_a(m2)");
      if ( what & PRT_D ) sptFprintf(F," depth(m)");
      if ( what & PRT_Z ) sptFprintf(F," surf_elev(m)");
      if ( what & PRT_XY) sptFprintf(F," xy-coordinates");
      sptFprintf(F,"\n");
    } else {
      sptFprintf(F,"*** id time(min)");
      if ( what & PRT_Q ) sptFprintf(F," flow(ft3/s)");
      if ( what & PRT_A ) sptFprintf(F," wet_a(ft2)");
      if ( what & PRT_D ) sptFprintf(F," depth(ft)");
      if ( what & PRT_Z ) sptFprintf(F," surf_elev(ft)");
      if ( what & PRT_XY) sptFprintf(F," xy-coordinates");
      sptFprintf(F,"\n");
    }
    // print the rest
    for (int jj=0; jj<_num_nodes; jj++) {
      if (!NAMES) sptFprintf(F, "%d %6.2f", _NS[jj]->Id(), tnow/60.0);
      else        sptFprintf(F, "%s %6.2f", NAMES[_NS[jj]->Id()], tnow/60.0);
      if ( what & PRT_Q) sptFprintf(F, PRT_FMT, _X[_NS[jj]->QIdx()] * fq);
      if ( what & PRT_A) sptFprintf(F, PRT_FMT, _X[_NS[jj]->AIdx()] * fa);
      if ( what & PRT_D) sptFprintf(F, PRT_FMT, _Depth[jj] * fd);
      if ( what & PRT_Z) sptFprintf(F, PRT_FMT, _Elevation[jj] * fd);
      if ( what & PRT_XY) sptFprintf(F, PRT_XY_FMT, _NS[jj]->X(), _NS[jj]->Y());
      sptFprintf(F,"\n");
    }
  } else {  // we only print at the frequency being asked for
    if ( tnowmin < pstart) return;


    if (  prev_t < 0 && pstart == 0 ) {
      // print the header
      if ( OPT.UseMetric() == 1 )  {
	sptFprintf(F,"*** id time(min)");
	if ( what & PRT_Q ) sptFprintf(F," flow(m3/s)");
	if ( what & PRT_A ) sptFprintf(F," wet_a(m2)");
	if ( what & PRT_D ) sptFprintf(F," depth(m)");
	if ( what & PRT_Z ) sptFprintf(F," surf_elev(m)");
	if ( what & PRT_XY) sptFprintf(F," xy-coordinates");
	sptFprintf(F,"\n");
      } else {
	sptFprintf(F,"*** id time(min)");
	if ( what & PRT_Q ) sptFprintf(F," flow(ft3/s)");
	if ( what & PRT_A ) sptFprintf(F," wet_a(ft2)");
	if ( what & PRT_D ) sptFprintf(F," depth(ft)");
	if ( what & PRT_Z ) sptFprintf(F," surf_elev(ft)");
	if ( what & PRT_XY) sptFprintf(F," xy-coordinates");

	sptFprintf(F,"\n");
      }

      for (int jj=0; jj<_num_nodes; jj++) {
	if (!NAMES) sptFprintf(F, "%d %d", _NS[jj]->Id(), 0 );
	else        sptFprintf(F, "%s %d", NAMES[_NS[jj]->Id()], 0);
	if ( what & PRT_Q) sptFprintf(F, PRT_FMT, _X[_NS[jj]->QIdx()] * fq);
	if ( what & PRT_A) sptFprintf(F, PRT_FMT, _X[_NS[jj]->AIdx()] * fa);
	if ( what & PRT_D) sptFprintf(F, PRT_FMT, _Depth[jj] * fd);
	if ( what & PRT_Z) sptFprintf(F, PRT_FMT, _Elevation[jj] * fd);
	if ( what & PRT_XY) sptFprintf(F, PRT_XY_FMT, _NS[jj]->X(), _NS[jj]->Y());
	sptFprintf(F,"\n");
      }
      // move print time stamp forward
      prev_t += interval;
    } else if ( tnowmin >= prev_t + interval )  {
      int tprt = prev_t + interval;
      while ( prev_t + interval <= tnowmin ) {
	if ( OPT.UseMetric() == 1 ) {
	  sptFprintf(F,"*** id time(min)");
	  if ( what & PRT_Q ) sptFprintf(F," flow(m3/s)");
	  if ( what & PRT_A ) sptFprintf(F," wet_a(m2)");
	  if ( what & PRT_D ) sptFprintf(F," depth(m)");
	  if ( what & PRT_Z ) sptFprintf(F," surf_elev(m)");
	  if ( what & PRT_XY) sptFprintf(F," xy-coordinates");
	  sptFprintf(F,"\n");
	} else {
	  sptFprintf(F,"*** id time(min)");
	  if ( what & PRT_Q ) sptFprintf(F," flow(ft3/s)");
	  if ( what & PRT_A ) sptFprintf(F," wet_a(ft2)");
	  if ( what & PRT_D ) sptFprintf(F," depth(ft)");
	  if ( what & PRT_Z ) sptFprintf(F," surf_elev(ft)");
	  if ( what & PRT_XY) sptFprintf(F," xy-coordinates");
	  sptFprintf(F,"\n");
	}
	for (int jj=0; jj<_num_nodes; jj++) {
	  if (!NAMES) sptFprintf(F, "%d %d", _NS[jj]->Id(), prev_t+interval);
	  else        sptFprintf(F, "%s %d", NAMES[_NS[jj]->Id()], prev_t+interval);
	  // in order for interpolation, we recompute depth and elevation at the previous
	  // time point on the fly. Note that "tprt" is given in minutes
	  if ( what & PRT_Q) sptFprintf(F, PRT_FMT, lin_interp(
							   _tm1, _Xtm1[_NS[jj]->QIdx()], 
							   _t,_X[_NS[jj]->QIdx()], tprt*60.0) * fq);
	  if ( what & PRT_A) sptFprintf(F, PRT_FMT, lin_interp(
							   _tm1,_Xtm1[_NS[jj]->AIdx()], 
							   _t,_X[_NS[jj]->AIdx()], tprt*60.0) * fa );
	  if ( what & PRT_D) sptFprintf(F, PRT_FMT, lin_interp(
							   _tm1, _NS[jj]->GetAbsoluteDepth(_Xtm1[_NS[jj]->AIdx()]), 
							   _t,_Depth[jj], tprt*60.0) * fd );
	  if ( what & PRT_Z) sptFprintf(F, PRT_FMT, lin_interp(
							   _tm1, _NS[jj]->GetElevation(_Xtm1[_NS[jj]->AIdx()]), 
							   _t,_Elevation[jj], tprt*60.0) * fd);

	  if ( what & PRT_XY) sptFprintf(F, PRT_XY_FMT, _NS[jj]->X(), _NS[jj]->Y());
	  sptFprintf(F,"\n");
	}
	// move print time stamp forward
	prev_t += interval;
	tprt += interval;
      }
    }


  }
}

/* 
   steady solve, this version uses forward Euler method to come to a close enough steady
   state before using the pseudo dynamic method
 
   it relies on the completed tchk results: Q at each node, and A at the outpouring node
   stored

   assuming Q is constant, we can compute the derivative dA/dx from steady-state SV. then
   use explicite Euler's scheme to compute A from downstream to upstream.

   current implementation assumes there is only ONE root node

*/
int Subcatchment::SteadySolve(int jac_num, int max_iter, double tol, int use_acc) {
  MyStack<int>      ST;
  MyStack<double>   ST_A;

  int       tmpbuf[16], status, cur, nxt, fo[32];
  Node      *pn;
  double    aa, qq, len, ratio[32];
  if ( _G.GetNumRoots() != 1 ) {
    fprintf(stdout,"[EE] Bummer: requires one (and only one) root node. Found %d\n",
	    _G.GetNumRoots() );
    return (-1);
  }
  
#if 0
  // option 1
  // DFS traversal + BE formula, deprecated

  double t1,t2, dadx, correction, naa;
  _G.GetRoots(&tmpbuf[0]); // too cautious?
  cur = tmpbuf[0];  

  pn = this->GetNode( cur );  

  aa = _G.GetA( cur );
  qq = _G.GetQ( cur );

  // we don't do anything for the root node, just copy the values
  pn->Q0() = qq;
  pn->A0() = aa;

  //  fprintf(stdout,"a=%.2f q=%.2f\n", aa, qq);

  // by checking the length, we would know whether we have a junction point or not
  status = _G.GetUpstreamLength( cur, len);
  if ( status != 0 ) {
    fprintf(stdout,"[EE] Bummer: the root node doesn't look right\n");
    return (-1);
  }
  
  nxt = _G.GetUpstreamNode( cur );
  // compute the values for the upstream nodes and push them into the stack
  t1 = OPT.Gravity() * (aa*pn->GetSR() - pn->Nsq()*qq*qq*pn->XS()->GetEqFriction(aa));
  t2 = OPT.Gravity() * aa * pn->XS()->GetDepthdA(aa) - qq*qq/(aa*aa);
  dadx = t1/t2;
  if ( dadx < 0 ) dadx=0.0;   // as a safe guard
  correction = 0.2*len*dadx;
  if ( aa > correction ) naa = aa - correction;
  else                   naa = aa;
  
  // get ready for the DFS traversal
  //   top of the stack:  ST: next 
  ST.Push( nxt );
  ST_A.Push( naa );

  // Standard DFS
  while (1) {
    status = ST.Pop(cur);
    if ( status == -1 ) break; 

    ST_A.Pop(aa);
    
    // store the values
    pn = this->GetNode( cur );
    qq = _G.GetQ( cur );
    pn->Q0() = qq;
    pn->A0() = aa;

#if 0
    fprintf(stdout,"a=%.3e q=%.3e\n", aa, qq);
#endif

    // compute for the upstream values and push to the stack
    status = _G.GetUpstreamLength( cur, len);
    if ( status == 0 ) {
      nxt = _G.GetUpstreamNode( cur );
      t1 = OPT.Gravity() * (aa*pn->GetSR() - pn->Nsq()*qq*qq*pn->XS()->GetEqFriction(aa));
      t2 = OPT.Gravity() * aa * pn->XS()->GetDepthdA(aa) - qq*qq/(aa*aa);
      dadx = t1/t2;
      if ( dadx < 0 ) dadx=0.0;
      correction = 0.2*len*dadx;
      
#if 0
      fprintf(stdout,"     correction=%.3e\n", correction);
#endif

      if ( aa > correction) naa = aa - correction;
      else                  naa = aa;

      ST.Push( nxt );
      ST_A.Push( naa );

    } else if ( status > 0 ) {
      _G.GetJunctionFanouts( cur, &fo[0], &ratio[0] );
      for (int jj=0; jj<status; jj++) {
	ST.Push( fo[jj] );
	ST_A.Push( aa * ratio[jj]);
      }
    } // else we are at leaves, simply do nothing
  }

#else
  // option 2
  // Use A associated with normal depth
  // no need to do DFS traversal, but use stack so that we can handle branches nicely

  _G.GetRoots(&tmpbuf[0]);
  cur = tmpbuf[0];  
  pn = this->GetNode( cur );  

  aa = _G.GetA( cur );
  qq = _G.GetQ( cur );

  double prev_aa = aa;

  // we don't do anything for the root node, since A is already given, just copy the values
  pn->Q0() = qq;
  pn->A0() = aa;

  nxt = _G.GetUpstreamNode( cur );
  ST.Push( nxt );
  
  // use stack
  while (1) {
    status = ST.Pop(cur);
    if ( status == -1 ) break; // stack is empty, done!
    
    pn = this->GetNode( cur );
    qq = _G.GetQ( cur );
    pn->Q0() = qq;

    aa = pn->KinematicEstimate( qq );
    if ( aa < 0 ) { // failed to find A, use a conservative value;
      pn->A0() = 5.0*prev_aa;
    } else {
      pn->A0() = aa;
      prev_aa = aa;
    }

#if 0
    fprintf(stdout,"a=%.3e q=%.3e\n", aa, qq);
#endif

    status = _G.GetUpstreamLength( cur, len); // whether we are at branch point or not
    if ( status == 0 ) {
      nxt = _G.GetUpstreamNode( cur );
      ST.Push( nxt );
    } else {
      _G.GetJunctionFanouts( cur, &fo[0], &ratio[0] );
      for (int jj=0; jj<status; jj++) {
	ST.Push( fo[jj] );
      }
    }
  }

#endif
  
  // copy to X
  this->InitSolutions();

  int rc=0;
  
#if 0
  // do Newton's iterations to solve for A only
  // by skipping it, the second phase (pseudo time marching will have to spend many more
  // iterations to achieve convergence
  
  if (use_acc==1) {
    int dobounding=1, num_iter=0;
    double maxnorm, minnorm, tdiff;
    FILE *F=NULL;

    // negative step size triggers solveing for dynamic equations only
    rc = NonLinearStep(0.0, -1.0, jac_num, max_iter, tol, dobounding, 
		       num_iter, minnorm, maxnorm, tdiff, F);

    if ( rc < 0 ) { // failed to converge
      fprintf(stdout,"Bummer: steady state solve failed to converge\n");
    } else {
      fprintf(stdout,"[II]: Steady state solve converged in %d steps\n", num_iter);
    }
  }
#endif
  return (rc);
}

// steady solve: this version uses pseudo dynamic method. It can be slow.
int Subcatchment::SteadySolve(double init_dt, int jac_num, int max_iter, double tol) {
  int rc, cnt, t_num_iter, num_iter;
  double tdif, dt, dt_used, pt, spin_r;
  int   negcounter = 0, poscounter=0;
  int   dc_iter=0;
  double maxnorm, minnorm;
  const int dobounding=1;
  FILE *F = NULL;

  cnt = 0;
  t_num_iter = 0;

  // spin-up for lateral flows
  if ( OPT.SpinUpTime() > 0 ) {
    dt = init_dt;
    pt = 0.0;
    while ( pt <= ((double) OPT.SpinUpTime()) ) {
      spin_r = pt/((double)OPT.SpinUpTime());
      if ( spin_r>1.0 ) spin_r = 1.0;
      if ( spin_r<=0.0) spin_r = 0.0;
      MasStvEqn::SetSpinCoef( spin_r );

      rc =  NonLinearStep(0.0, dt, jac_num, max_iter, tol, dobounding, 
			  num_iter, minnorm, maxnorm, tdif, F);
      t_num_iter += num_iter;
      dc_iter++;

      dt_used = dt;
      if ( rc < 0 ) {  // failed to converge

	STAT.Add_ss_step(0.0, dt, num_iter, -1);
	
	dt *= 0.5;
	negcounter++;
	if ( dt < OPT.MinDT() ) {
	  fprintf(stdout,"Bummer, Steady State solve failed (dt=%.2e)\n", dt);
	  cnt = -1;
	  return cnt;
	}
	poscounter = 0;

	if (F) {
	  fprintf(F, "convg=0;\n");
	  fclose(F);
	}
	
      } else {   // converged, determine the new step size

	pt += dt_used;

	STAT.Add_ss_step(0.0, dt, num_iter, 1);
	
	if ( num_iter < 8 ) poscounter++;
	else                poscounter=0;
	
	if (poscounter > 5 ) {
	  dt_used = 2*dt;
	  poscounter = 0;
	  if ( dt_used < 50.0 ) dt = dt_used;
	  else dt = 50.0;
	}
	
	if (F) {
	  fprintf(F, "convg=1;\n");
	  fclose(F);
	}
	
      }
      
      if (OPT.DebugLevel() >= 1 ) {
	CalSolStat();
	fprintf(stdout,"  S Iter w/ dt=%.2e %s, nxt_dt=%.2e, num_iter=%3d, max_norm=%.2e, tdif=%.2e, n_sc=%2d, max_f=%.3f\n",
		dt_used, rc<0?"FAILED":"passed", dt, num_iter, maxnorm, tdif, GetNumSuperC(), GetMaxFroude());
      }
      
      cnt++;
    }
    fprintf(stdout, "[II]: Spin-up completed for %d seconds, with %d steps and %d iterations.\n", OPT.SpinUpTime(), cnt, t_num_iter);
  }

  // reset
  dt = init_dt;
  tdif = OPT.BogusLimit();
  MasStvEqn::SetSpinCoef( 1.0 );
  
  while ( tdif > OPT.SSTol() ) {
    
#ifdef DBGSTEADY
    char fn[512];
    char vn[512];
    sprintf(fn,"dc_%d.m",dc_iter);
    F = fopen(fn,"w");
#endif

    rc =  NonLinearStep(0.0, dt, jac_num, max_iter, tol, dobounding, 
			num_iter, minnorm, maxnorm, tdif, F);
    t_num_iter += num_iter;
    dc_iter++;

    dt_used = dt;
    if ( rc < 0 ) {  // failed to converge

      STAT.Add_ss_step(0.0, dt, num_iter, -1);

      dt *= 0.5;
      negcounter++;
      if ( dt < OPT.MinDT() ) {
	fprintf(stdout,"Bummer, Steady State solve failed (dt=%.2e)\n", dt);
	cnt = -1;
	return cnt;
      }
      poscounter = 0;

      if (F) {
	fprintf(F, "convg=0;\n");
	fclose(F);
      }

    } else {   // converged, determine the new step size

      STAT.Add_ss_step(0.0, dt, num_iter, 1);

      if ( num_iter < 8 ) poscounter++;
      else                poscounter=0;

      if (poscounter > 5 ) {
	dt_used = 2*dt;
	poscounter = 0;
	if ( dt_used < 50.0 ) dt = dt_used;
	else dt = 50.0;
      }

      if (F) {
	fprintf(F, "convg=1;\n");
	fclose(F);
      }

    }

    if (OPT.DebugLevel() >= 1 ) {
      CalSolStat();
      fprintf(stdout,"  S Iter w/ dt=%.2e %s, nxt_dt=%.2e, num_iter=%3d, max_norm=%.2e, tdif=%.2e, n_sc=%2d, max_f=%.3f\n",
	      dt_used, rc<0?"FAILED":"passed", dt, num_iter, maxnorm, tdif, GetNumSuperC(), GetMaxFroude());
    }

    cnt++;
  }

  fprintf(stdout, "[II]: Number of nodes = %d, size of matrix = %d\n", _num_nodes, _num_eqns);
  fprintf(stdout, "[II]: SteadyState phase II converged in %d steps with %d iterations.\n", cnt, t_num_iter);

  return cnt;
}

// query to find the closet break point
double Subcatchment::QueryBreakPoint(double tnow, int &atbrkpt) {
  int k;
  double nt, t;
  QrSrcEqn *Q;
  int rc=0;

  atbrkpt = 0;
  t = tnow + FLT_MAX;
  for (k=0; k<_num_eqns;k++) {
    if (_EQ[k]->Type() != Q_SRC ) continue;
    Q = (QrSrcEqn*)_EQ[k];
    rc += Q->Src()->GetNextT(tnow, nt);
    //    nt = Q->Src()->NextT(oldt);
    t = MIN(t, nt);
  }
  atbrkpt = rc>0 ? 1 : 0;
 
  return (t-tnow);
}

/* unsteady solve, returns the number of steps if succeeded, -1 if failed
   if successful and OUT!=NULL, print the values to OUT at each time point
 */
int Subcatchment::UnsteadySolve(double final_t, int jac_num, int max_iter, double tol, 
				Waveforms *WV, sptFile OUT, int what2print, int pstart, char **NAMES) {

  double current_t = _t;  // we can call this many times, pick up t from previous runs
  int    rc;
  double optdt, dtused;
  int    t_num_iter, num_iter,cnt;
  int    niter_here;
  int    first_start;
  int    poscounter, negcounter;
  double maxnorm, minnorm, tdif;
  double tf_hr = final_t/3600;
  const int do_bounding=0;         // we don't need bounding for unsteady
  
  cnt = 0;
  optdt = 200.0;
  dtused = 15.0;    // initial value 15 seconds
  t_num_iter = 0;
  poscounter = 0;
  negcounter = 0;

  if ( final_t < _t ) return 0; // no need to do anything 

  first_start = ( _dt > 0 );

  if (WV) WV->Init(_num_nodes);

  // restore Xp from the previous time point
  memcpy(&_Xp[0], &_Xtm1[0], _num_eqns*sizeof(double));

  // print a header
  if (OUT!=NULL) {

    sptFprintf(OUT,"*** SPRNt Results. Netlist: \"%s\" Epoch: \"%s\" Reqested Tstop: ", STAT.InFile(), STAT.Epoch() );

    if ( (((int)OPT.StopTime())/3600) > 0 ) {
      sptFprintf(OUT,"%d hr %d min\n", ((int)OPT.StopTime())/3600, (((int)OPT.StopTime())/60)%60 );
    }  else  {
      sptFprintf(OUT,"%d min\n",((int)OPT.StopTime())/60);
    }

  }

  while ( 1 ) {
    
    /* determine the time step we should use */
    optdt = CalFroude();         // optimal time step based on Froude #

    if ( first_start ) { // case1: first time enter restart, use previous 
      dtused = _dt;
      first_start = 0;
    } else if ( OPT.FixedStep() > 0 ) { // case2: fixed time step
      dtused = OPT.FixedStep(); 
    } else {                     // otherwise we tune the time step
      if ( poscounter > 4 ) {    // not too many iterations, double
	dtused *= 2.0;
	poscounter = 0;
      }
      if ( dtused > optdt ) dtused = optdt;  // but we don't go too far
    }

    if (OPT.DebugLevel() >= 2 ) {
      fprintf(stdout, "    next proposed dt=%.2e\n", dtused);
    }

    rc = -1;
    negcounter = 0;
    niter_here = 0;

    // nonlinear iteration, 
    while ( rc < 0 ) {
      rc =  NonLinearStep(current_t, dtused, jac_num, max_iter, tol,do_bounding,
			  num_iter, minnorm, maxnorm, tdif);
      t_num_iter += num_iter;
      niter_here += num_iter;
      if ( rc < 0 ) {

	STAT.Add_dy_step(current_t, dtused, num_iter, -1);

	dtused *= 0.5;
	negcounter++;
	if ( OPT.DebugLevel() >= 2) {
	  fprintf(stdout,"   reduce time step to %.2e\n", dtused);
	}
	if ( dtused < OPT.MinDT() ) {
	  fprintf(stderr, "[EE] Bummer: time step too small at time point %.4e\n", current_t);
	  return ERROR;
	}
      }
    }

    // NL has converged

    STAT.Add_dy_step(current_t, dtused, num_iter, -1);

    // we tried to maintain the number of iterations between 14 and 25
    if ( negcounter == 0 && niter_here < 14 ) poscounter++; // magic number
    else poscounter = 0;

    if (OPT.DebugLevel() >= 1 ) {
      int nbk = CheckBankFull();
      int nmn = CheckMinimalA();
      fprintf(stdout,"  U time=%.4e (%6.2f out of %5.2f hrs), num_iter=%2d, falter=%d, dt=%.2e, n_sc=%3d, max_f=%5.2f, n_bk=%3d, n_mn=%3d\n", 
	      _t, _t/3600, tf_hr, niter_here, negcounter, dtused, GetNumSuperC(), GetMaxFroude(), nbk, nmn );
    }

    // compute depth and elevation when needed
    if ( (what2print & PRT_D) || (what2print & PRT_Z) || WV ) ComputeDepthAndElevation(what2print);

    if ( WV ) {
      int rc= WV->Advance();
      if ( rc < 0 ) {
	fprintf(stdout,"[WW]: Error allocating waveform storage, no more results stored.\n");
      } else {
	WV->StoreQplusA(current_t, _X);
	WV->StoreD(current_t, &_Depth[0]);
	WV->StoreZ(current_t, &_Elevation[0]);
      }
    }

    if (OUT!=NULL ) {
      PrintOnDemand(OUT, _t, pstart, what2print, NAMES);
    }

    current_t += dtused;  // next time point
    _dt = dtused;         // record the past used dt

    if ( _t >= final_t ) {
      break;
    }

    cnt++;
  }

  fprintf(stdout, "[II]: UnsteadyState finished w/ %d time points, %d iterations.\n", cnt, t_num_iter);

  return cnt;
}

// returns 0 if converged
// return -1 if failed
int Subcatchment::NonLinearStep(double t, double dt, int jac_num, int max_iter, double tol, int do_bounding,
				int &num_iter, double &minn, double &maxn, double &tdiff, FILE *F) {
  int rc = -1;
  int cnt;
  double l2norm, dummy1, dummy2;
  double neg_alpha = -OPT.Alpha();
  const int onei = 1;
  const int twoi = 2;

  assert ( jac_num > 0 );  // otherwise we are going dodo

  tdiff = OPT.BogusLimit();
  minn = 1e12;
  maxn = -1e12;

  // copy Xp to X as the initial value
  memcpy(&_X[0], &_Xp[0], _num_eqns*sizeof(double));


  for (cnt=0; cnt < max_iter; cnt++) {

    if (cnt < jac_num)  Evaluate(t, dt);
    else                EvaluateRHS(t, dt);


    l2norm = this->GetNorm();

    if (OPT.DebugLevel() >= 2 ) {
      fprintf(stdout,"    NL step %2d out of %2d, jac=%d, bnd=%1d, norm=%.2e\n", 
	      cnt, max_iter, jac_num, do_bounding, l2norm);
    }

    minn = MIN(minn, l2norm);
    maxn = MAX(maxn, l2norm);
    if ( l2norm < tol ) { 
      rc = 0;
      break;
    }

    // solve
    _s->Solve((double*)_RHS, (double*)_DX);

    if (do_bounding) {
      memcpy(&_Xhigh[0], &_X[0], _num_eqns*sizeof(double));
      memcpy(&_Xlow[0],  &_X[0], _num_eqns*sizeof(double));
    }

    FORTRAN(daxpy)(&_num_eqns, &neg_alpha, &_DX[0], &onei, &_X[0], &onei);

    if (do_bounding) {
      double lowlim = OPT.DownLimit();
      double highlim = OPT.UpLimit();
      FORTRAN(dscal)(&_num_nodes, &lowlim,  &_Xlow[1],  &twoi); // we bound A, hence start
								// from 1, and skip every
								// other entry
      FORTRAN(dscal)(&_num_nodes, &highlim, &_Xhigh[1], &twoi);
      for (int j=1; j<_num_eqns; j+=2) _X[j] = _X[j] < _Xlow[j]  ? _Xlow[j]  : _X[j];
      for (int j=1; j<_num_eqns; j+=2) _X[j] = _X[j] > _Xhigh[j] ? _Xhigh[j] : _X[j];
    }  

#ifdef DBGNL
    if (F) {
      char vname[512];
      sprintf(vname,"jac%d",cnt);
      _M.MatlabDump(F,vname);
      sprintf(vname,"X%d",cnt);
      MatlabDumpX(F,vname);
      sprintf(vname,"DX%d",cnt);
      MatlabDumpDX(F,vname);
      sprintf(vname,"RHS%d",cnt);
      MatlabDumpRHS(F,vname);
    }
#endif

    if ( CheckMinA(OPT.MinA()) == -1) break;
  }

  if ( rc == 0 ) tdiff = Accept(t, 1, dummy1, dummy2);

  num_iter = cnt;
  return rc;
}

int Subcatchment::CheckBankFull() {
  int  rc = 0;
  for (int j=0; j<_num_nodes; j++) rc += _NS[j]->XS()->ReachedBankFull( _X[ _NS[j]->AIdx() ] );    
  return rc;
}
int Subcatchment::CheckMinimalA() {
  int rc = 0;
  for (int j=0; j<_num_nodes; j++) rc += _NS[j]->XS()->ReachedMinimalA( _X[ _NS[j]->AIdx() ] );    
  return rc;
}

int Subcatchment::ComputeDepthAndElevation(int flags) {
  if ( flags & PRT_D ) {
    for (int jj=0; jj<_num_nodes; jj++) _Depth[jj] = _NS[jj]->GetAbsoluteDepth(_X[_NS[jj]->AIdx()]);
  }

  if ( flags & PRT_Z ) {
    for (int jj=0; jj<_num_nodes; jj++) _Elevation[jj] = _NS[jj]->GetElevation(_X[_NS[jj]->AIdx()]);
  }

  return 0;
}

// End
