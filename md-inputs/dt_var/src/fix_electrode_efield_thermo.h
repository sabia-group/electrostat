/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (GU), Kamila Savvidi (TUHH), Robert Meissner (Hereon, TUHH)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

// clang-format off
FixStyle(electrode/efield/thermo, FixElectrodeEfieldThermo);
// clang-format on

#else

#ifndef LMP_FIX_ELECTRODE_EFIELD_THERMO_H
#define LMP_FIX_ELECTRODE_EFIELD_THERMO_H

#include "fix_electrode_efield.h"

namespace LAMMPS_NS {
// forward decls

class ElectrodeTaglist;

class FixElectrodeEfieldThermo : public FixElectrodeEfield {

 public:
  FixElectrodeEfieldThermo(class LAMMPS *, int, char **);
  ~FixElectrodeEfieldThermo() override;
  double compute_vector(int) override;
  void post_constructor() override;

 protected:
  virtual void set_efield() override;

 private:
  class RanMars *thermo_random;
  double efield;
  double de_pol, de_measure;
  int efield_step;
};

}    // namespace LAMMPS_NS

#endif
#endif
