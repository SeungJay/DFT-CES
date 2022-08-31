/* -*- c++ -*- ----------------------------------------------------------
   DFT-CES core subroutines. Written by H.-K. Lim
   Copyright (C) 2016 M-design group @ KAIST
------------------------------------------------------------------------- */

#ifndef LMP_GRID_H
#define LMP_GRID_H

#include "pointers.h"

namespace LAMMPS_NS {

class Grid : protected Pointers {
 public:
  int gnx, gny, gnz;                   // # of grid points in 3 dim
  double gx[3],gy[3],gz[3];            // grid spacing vectors of (unit: Ang)
  int natoms;                          // # of atoms
  int *atomns;                         // atomic numbers for each atom
  double **basis;                      // cartesian coords of each atom (unit: Ang)
                                       // within unit cell (0 <= coord < 1)
  double *gvin;                        // grid values from QM pot (unit: Ry)
  double *gvout, *gvout_all;           // grid values for MD rho  (unit: Ry)

  Grid(class LAMMPS *, int, char **);
  ~Grid();

  void read_header(char *);
  void read_content(char *);
  void save_grid(char *, int);
};

}

#endif

