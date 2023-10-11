/* -*- c++ -*- ----------------------------------------------------------
   DFT-CES core subroutines. Written by H.-K. Lim
   Copyright (C) 2016 M-design group @ KAIST
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(gridforce,FixGridForce)

#else

#ifndef LMP_FIX_GRIDFORCE_H
#define LMP_FIX_GRIDFORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGridForce : public Fix {
 public:
  FixGridForce(class LAMMPS *, int, char **);
  ~FixGridForce();
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  double compute_scalar();
  void triInter(double, double, double, double, double*);
  double triInterValue(double,double,double,double,double,double,double,double,double,double,double);
  double gridValue(double*, int, int, int);
  void gridGrad(double*, int, int, int, double*);

 private:
  double weight;      // weight factor for grid energy and force, usually -1 for QE coupling
  int sfactor;        // supercell factor
  double tmp3[3], tmp4[4];
 
  int force_flag;

  double Egrid, Egrid_all;

};

}

#endif
#endif
