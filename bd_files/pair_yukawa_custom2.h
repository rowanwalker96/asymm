/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(yukawacustom2,PairYukawaCustom2);
// clang-format on
#else

#ifndef LMP_PAIR_YUKAWA_CUSTOM2_H
#define LMP_PAIR_YUKAWA_CUSTOM2_H

#include "pair.h"

namespace LAMMPS_NS {

class PairYukawaCustom2 : public Pair {
 public:
  PairYukawaCustom2(class LAMMPS *);
  virtual ~PairYukawaCustom2();
  virtual void compute(int, int);
  void init_style();
  void settings(int, char **);
  void coeff(int, char **);
  virtual double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  void *extract(const char *, int &);
  virtual double single(int, int, int, int, double, double, double, double &);

 protected:
  double cut_global;
  double kappa1;
  double kappa2;
  double hamaker;
  double *rad;
  double **cut, **a, **b, **offset;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

*/
