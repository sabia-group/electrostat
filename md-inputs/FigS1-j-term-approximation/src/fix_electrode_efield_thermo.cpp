/* ----------------------------------------------------------------------
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

#include "fix_electrode_efield_thermo.h"

#include "atom.h"
#include "charge_solver.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "math_const.h"
#include "modify.h"
#include "neigh_list.h"
#include "pair.h"
#include "random_mars.h"
#include "update.h"
#include "variable.h"

#include <cassert>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

FixElectrodeEfieldThermo::FixElectrodeEfieldThermo(LAMMPS *lmp, int narg, char **arg) :
    FixElectrodeEfield(lmp, narg, arg)
{
  vector_flag = 1;
  size_vector = 3;
  efield_step = -1;
  if (thermo_seed < 0) error->all(FLERR, "temp keyword not set for {}", style);
  thermo_random = new RanMars(lmp, thermo_seed);
  de_pol = de_measure = 0.;
}

/* ---------------------------------------------------------------------- */

double FixElectrodeEfieldThermo::compute_vector(int i)
{
  assert(efield_step == update->ntimestep);
  switch (i) {
    case 0:
      return efield;
    case 1:
      return de_pol;
    case 2:
      return de_measure;
    default:
      error->all(FLERR, "vector index {} not implemented for {}", i + 1, style);
  }
}

/* ----------------------------------------------------------------------
    Setup a fix efield for forces on the electrolyte
------------------------------------------------------------------------- */

void FixElectrodeEfieldThermo::post_constructor()
{
  assert(field.style != VarStyle::UNSET);
  // check for other efields and warn if found
  if ((modify->get_fix_by_style("^efield").size() > 0) && (comm->me == 0))
    error->warning(FLERR, "Other efield fixes found -- please make sure this is intended!");

  std::string var_efield = fixname + "_efield";
  input->variable->set(fmt::format("{} equal f_{}[1]", var_efield, fixname));
  std::string efield_call = fixname + "_efield all efield 0.0 0.0 v_" + var_efield;
  modify->add_fix(efield_call, 1);    // call fix command
}

/* ----------------------------------------------------------------------
    Update efield according to electrostat equation
------------------------------------------------------------------------- */

void FixElectrodeEfieldThermo::set_efield()
{
  assert(field.style != VarStyle::UNSET);
  assert(!symm_split);    // TODO remove symm_split keyword
  const double efield_0 = field.style == VarStyle::CONST ? field.value
                                                         : input->variable->compute_equal(field.id);
  const double vac_cap = charge_solver->vacuum_capacitance();
  if (update->ntimestep != update->beginstep) {    // TODO set efield if not starting with beginstep
    const double dt = update->dt;
    const double exp_dt = exp(-dt / thermo_time);
    const double noise_factor =
        1. / domain->prd[2] * sqrt((thermo_temp / vac_cap) * (1. - exp_dt * exp_dt));

    // polarization term
    double current = 0.;
    for (int i = 0; i < atom->nlocal; i++) {
      if (!(atom->mask[i] & groupbit)) current += atom->q[i] * atom->v[i][2];
    }
    MPI_Allreduce(MPI_IN_PLACE, &current, 1, MPI_DOUBLE, MPI_SUM, world);
    const double epsilon0 = evscale / MY_4PI;
    const double vol = domain->xprd * domain->yprd * domain->zprd;
    de_pol = -thermo_time * (1 - exp_dt) * current / (epsilon0 * vol);

    double de_decay = (efield - efield_0) * exp_dt;
    double de_random = thermo_random->gaussian() * noise_factor;
    efield = efield_0 + de_decay + de_pol + de_random;

    de_measure = (charge_solver->get_field() - efield_0) * exp_dt - de_decay;
  } else {
    efield = charge_solver->get_field();
  }
  charge_solver->set_efield(efield);
  efield_step = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

FixElectrodeEfieldThermo::~FixElectrodeEfieldThermo()
{
  delete thermo_random;
}
