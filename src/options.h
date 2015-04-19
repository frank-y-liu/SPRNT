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

/**************************************************************************
   options and constants
**************************************************************************/

#ifndef _OPTIONS_H
#define _OPTIONS_H

#include "flexvec.h"

class Options {
 private:
  double    _g;          // gravity
  double    _sqrt_g;     // sqrt of g
  double    _ftom;       // feet to meters
  double    _f3tom3;     // cubic feet to cubic meters
  double    _f2tom2;     // square feet to square meters
  int       _use_metric; // use metric or US
  double    _epsilon;    // for checking
  double    _beta;       // for lateral correction
  double    _minq;       // min Q
  double    _mina;       // minimal A
  double    _alpha;      // damping coefficient
  double    _max_c;      // max celerity
  double    _super_c;    // where to define super critical
  double    _tol;        // tolerance for Newton
  double    _ss_tol;     // tolerance for steady-state, relative
  double    _k_ratio;    // factor to emphasis ratio between internal and connection nodes
  double    _min_dt;     // minimal dt when we declare failure
  double    _dt_factor;  // factor to scale down dt when convergence is an issue
  double    _cele_r;     // max dt ratio
  double    _down_limit; // down limit for bounding
  double    _up_limit;   // up limit for bouding
  double    _bogus_limit;// bogus limit
  double    _ht;         // limit to envoke n correction, in meters, set to zero to
			 // disable it
  int       _check_only; // flag to run kinematic checking only
  int       _steady_acc; // flag to accelerate steady solve
  int       _spin_up;    // source ramping time (in seconds), set to <=0 to disable it
  int       _print_int;  // interval to print, in minutes, if 0, print everything
  int       _print_start; // starting time point to start printing, in minutes.
  int       _debug_level;// debug level, controls how much blabhing will be printed
  double    _epsilona;   // algorithmic tuning parameter for A
  double    _minb;       // value to guard the bottom so that we don't have P overflow
  double    _fixed_step; // fixed time step, if zero, use automatic time step
  double    _stop_time;  // stopping time
  int       _print_q;    // flag to print Q
  int       _print_a;    // flag to print A
  int       _print_d;    // flag to print D
  int       _print_z;    // flag to print surf elev
  int       _print_xy;   // flag to print xy coordinates
  double    _lmax;       // L max, for checking only
  double    _lmin;       // L min for checking only
  double    _min_n;      // minimal manning's N, for checking only
  char      _buf[512];    // for others to copy the content to
  // char      _ssfile[512]; // steady state file
  // char      _chksum[10];  // stores the checksum for both node names and Q/A/lateral
  // 			  // sources
  // char      _epoch[128];  // epoch for time 0

  void _defaults() {
    _g = 9.80665;
    _sqrt_g = 3.13156;
    _ftom = 0.3048;
    _f3tom3 = 0.028317;
    _f2tom2 = 0.092903;
    _use_metric = 0;
    _epsilon = 1e-6;
    _beta = 1.03;
    _minq = 1.0e-5;
    _mina = 1.0e-9;
    _alpha = 0.86;
    _max_c = 0.90;
    _super_c = 0.99;
    _tol = 1e-6;
    _ss_tol = 1e-8;  
    _k_ratio = 1e0;
    _min_dt = 1e-3;
    _dt_factor = 0.5;
    _cele_r = 85;
    _down_limit = 0.5;
    _up_limit = 4.0;
    _bogus_limit = 1e9;
    _ht = 0.08;
    _debug_level = 0;
    _epsilona = 1e-4;
    _minb = 0.1;
    _fixed_step = 0;
    _print_int = 0;
    _print_start = 0;
    _stop_time = 0.0;
    _print_q = 0;
    _print_a = 0;
    _print_d = 0;
    _print_z = 0;
    _print_xy = 0;
    _lmax = 9.9e+9;
    _lmin = 0.0;
    _min_n = 1e-3;
    _check_only = 0;
    _steady_acc = 1;
    _spin_up = 0;
    _buf[0]= 0;
    // _ssfile[0]=0; // turn off by default
    // _chksum[0] = 0;
    // strcpy(_epoch,"1970-01-01T00:00:00Z");
  }

 public:
  Options() { _defaults();}
  ~Options() {}

  /* constants */
  double Gravity() const { return _g; }
  double SqrtG() const { return _sqrt_g; }
  double FtoM() const { return _ftom; }
  double F2toM2() const { return _f2tom2; }
  double F3toM3() const { return _f3tom3; }

  /* coefficients, can be changed */
  int& UseMetric() { return _use_metric; }
  double& Epsilon() { return _epsilon; }
  double& Beta() { return _beta; }
  double& MinQ() { return _minq; }
  double& MinA() { return _mina; }
  double& Alpha() { return _alpha;}
  double& MaxC() { return _max_c;}
  double& SuperC() { return _super_c; }
  double& Tol() { return _tol;}
  double& SSTol() { return _ss_tol; }
  double& KRatio() { return _k_ratio;}
  double& MinDT()  { return _min_dt;}
  double& DtFactor() { return _dt_factor;}
  double& MaxDtR() { return _cele_r; }
  double& DownLimit() { return _down_limit; }
  double& UpLimit() { return _up_limit; }
  double BogusLimit() const { return _bogus_limit; }  // this is read only
  double& HT() { return _ht; }
  int& DebugLevel() { return _debug_level; }
  double &EpsilonA() { return _epsilona; }
  double &MinB() { return _minb; }
  double &FixedStep() { return _fixed_step; }
  double &StopTime() { return _stop_time; }
  int& PrintInterval() { return _print_int; }
  int& PrintStart() { return _print_start; }
  int& CheckOnly() { return _check_only; }
  int& SteadyAcc() { return _steady_acc; }
  int& SpinUpTime() { return _spin_up; }

  int& PrintQ() { return _print_q; }
  int& PrintA() { return _print_a; }
  int& PrintD() { return _print_d; }
  int& PrintZ() { return _print_z; }
  int& PrintXY() {return _print_xy;}
  double& LMax() { return _lmax; }
  double& LMin() { return _lmin; }
  double& MinN() { return _min_n; }
  char* BufStr() { return _buf; }
  void  CopyToBuf(char *s) { assert(s!=NULL); strncpy(_buf, s,512); _buf[511]=0;}
  
  // char* SSFile() { if (_ssfile[0] != 0 ) return _ssfile; else return NULL; }
  // void  SetSSFile(char *s) { assert(s!=NULL); strncpy(_ssfile,s,512); _ssfile[511]=0;}

  // char* ChkSum() { if (_chksum[0]!=0) return _chksum; else return NULL; }
  // void  SetChkSum(char *s) { assert(s!=NULL); strncpy(_chksum, s, 8); _chksum[8]=0;}
  // char* Epoch()  { return _epoch; }
  // void  SetEpoch(const char *s) { assert(s!=NULL); strncpy(_epoch, s, 128); _epoch[127]=0; }
};

class Stats {
  typedef struct {
    double     _time;
    double     _step;
    int        _num_iter;  // negative means not converged, positive means converged
  } Behavior;
  
 private:
  long int      _n_nodes;
  long int      _n_segs;
  long int      _n_juncs;
  long int      _n_qsrcs;
  long int      _n_lsrcs;
  long int      _n_bnds;

  long int      _ss_steps;
  long int      _dy_steps;

  GrowVec<Behavior, 2048>   _steady_stats;
  GrowVec<Behavior, 8192>   _dynamic_stats;

  // stats related to a particular netlist
  char           _ssfile[512];
  char           _in_file[512];
  char           _chksum[10];
  char           _epoch[128];

 public:
 Stats():_n_nodes(0),_n_segs(0),_n_juncs(0),_n_qsrcs(0),_n_lsrcs(0),_n_bnds(0) {
    _ss_steps = 0;
    _dy_steps = 0;
    _ssfile[0] = 0;  // turned off
    _in_file[0] = 0;
    _chksum[0] = 0;
    strcpy(_epoch, "1970-01-01T00:00:00Z");
  }
  
  ~Stats() {}

  // methods
  long int& N_Nodes() { return _n_nodes; }
  long int& N_Segs() { return _n_segs; }
  long int& N_Juncs() { return _n_juncs; }
  long int& N_Qsrcs() { return _n_qsrcs; }
  long int& N_Lsrcs() { return _n_lsrcs; }
  long int& N_Bnds() { return _n_bnds; }

  const long int N_ss_steps() const { return _ss_steps; }
  const long int N_dy_steps() const { return _dy_steps; }
  
  void Add_ss_step(double t, double step, int niter, int isconverged) {
    int sg = isconverged ? 1 : -1;
    _steady_stats[(int)_ss_steps]._time = t;    // down converting, could be dangerous
    _steady_stats[(int)_ss_steps]._step = step;
    _steady_stats[(int)_ss_steps]._num_iter = sg*niter;
    _ss_steps++;
  }

  void Add_dy_step(double t, double step, int niter, int isconverged) {
    int sg = isconverged ? 1 : -1;
    _dynamic_stats[(int)_dy_steps]._time = t;
    _dynamic_stats[(int)_dy_steps]._step = step;
    _dynamic_stats[(int)_dy_steps]._num_iter = sg*niter;
    _dy_steps++;
  }

  char* SSFile() { if (_ssfile[0] != 0 ) return _ssfile; else return NULL; }
  void  SetSSFile(char *s) { assert(s!=NULL); strncpy(_ssfile,s,512); _ssfile[511]=0;}
  
  char* ChkSum() { if (_chksum[0]!=0) return _chksum; else return NULL; }
  void  SetChkSum(char *s) { assert(s!=NULL); strncpy(_chksum, s, 8); _chksum[8]=0;}
  
  char* InFile() { if (_in_file[0]!=0) return _in_file; else return NULL; }
  void SetInFile(const char *s) { assert(s!=NULL); strncpy(_in_file,s,128); _in_file[127]=0;}

  char* Epoch()  { return _epoch; }
  void  SetEpoch(const char *s) { assert(s!=NULL); strncpy(_epoch, s, 128); _epoch[127]=0; }

};

extern Options OPT;
extern Stats   STAT;

#endif

// Local Variables:
// mode: c++
// End:

