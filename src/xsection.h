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
  cross section definitions
 ********************************************************************/

#ifndef _XSECTION_H
#define _XSECTION_H

#include <math.h>
#include "base.h"
#include "flexvec.h"

/* virtualized so that we can have multiple implementations */
class XSection {
 protected:
  int          _id;
  float        _aa_min;
  float        _aa_bf;

 public:
 XSection(int id):_id(id),_aa_min(1e-6),_aa_bf(9.99e6) {} // arbitrary
  virtual ~XSection() {}

  // methods
  int GetId() const { return _id; };
  virtual float GetMinA() const{ return _aa_min;}
  virtual float GetBkfA()  const { return _aa_bf;}

  virtual double GetDepth(double a) {return 0.0;}
  virtual double GetDepthdA(double a) { return 0.0;}
  virtual double GetWidth(double a) {return 0.0;}
  virtual double GetCentroid(double a) {return 0.0;}
  virtual double GetCentroiddA(double a) {return 0.0;}
  virtual void   GetCentroid(double a, double& cen, double& cenda) {} //they might be
								      //cheaper
  // EqFriction term is basically P^(4/3)/ A^(7/3)
  virtual double GetEqFriction(double a) {return 0.0;}
  virtual double GetEqFrictiondA(double a) {return 0.0;}
  virtual void   GetEqFriction(double a, double &ef, double &efda) {}
  virtual void   GetHydroRadius(double a, double &r, double &drda) {}

  // finite difference method to check derivatives, they are not intended for real code
  virtual double CheckCentroiddA(double a) { return 0.0;}
  virtual double CheckEqFrictiondA(double a) { return 0.0;}

  virtual double GetAbyDepth(double h) { return 0.0; }

  virtual int ReachedBankFull(double a) { return 0; }
  virtual int ReachedMinimalA(double a)     { return 0; }

};

/****************************************************
  rectangular cross section
****************************************************/
class R_XSection : public XSection {   
 private: 
  double      _b0;    // width
  
 public:
   R_XSection(int id, double b0):XSection(id) {assert(b0>0); _b0=b0;}
  ~R_XSection() {};

  double GetDepth(double a);
  double GetDepthdA(double a);
  double GetWidth(double a);
  double GetCentroid(double a);
  double GetCentroiddA(double a);
  double GetEqFriction(double a);
  double GetEqFrictiondA(double a);

  void GetCentroid(double a, double &cen, double &cenda);
  void GetEqFriction(double a, double &ef, double &efda);
  void GetHydroRadius(double a, double &r, double &drda);

  double GetB0() const { return _b0; }

  double GetAbyDepth(double h) { return (_b0*h); }

};

/****************************************************
 trapezoidal cross section
****************************************************/
class T_XSection : public XSection {
 private:
  double       _b0;     // bottom width
  double       _effb0;  // scaled bottom width, b0/s/2
  double       _s;      // sidewall slope

 public:
 T_XSection(int id):XSection(id) { _b0 = -1; _s=-1; }
 T_XSection(int id, double b0, double s):XSection(id){
    assert(b0>0 && s>0);
    _b0 = b0;
    _s = s;
    _effb0 = 0.5*b0/s;
  };
  ~T_XSection() {};

  void Setb0(double b) {assert(b>0); _b0=b;}
  void SetS(double s) {assert(s>0); _s=s;}

  double GetDepth(double a);
  double GetDepthdA(double a);
  double GetWidth(double a);
  double GetCentroid(double a);
  double GetCentroiddA(double a);
  double GetEqFriction(double a);
  double GetEqFrictiondA(double a);

  void GetCentroid(double a, double &cen, double &cenda);
  void GetEqFriction(double a, double &ef, double &efda);
  void GetHydroRadius(double a, double &r, double &drda);

  double GetAbyDepth(double h) { return (h*( h*_s + _b0) ); }

  double GetB0() const { return _b0; }
  double GetS() const { return _s; }
};

// XY cross section is defined in a separate file

#endif

// Local Variables:
// mode: c++
// End:
