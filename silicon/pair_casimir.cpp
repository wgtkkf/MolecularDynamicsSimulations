/* ----------------------------------------------------------------------
   pair_casimir.cpp   --  custom dispersion / Casimir pair style
------------------------------------------------------------------------- */

#include "pair_casimir.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairCasimir::PairCasimir(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  restartinfo = 0;     // set to 1 only if you also implement write_restart()
  single_enable = 1;
}

/* ---------------------------------------------------------------------- */

PairCasimir::~PairCasimir()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(c6);
    memory->destroy(lam);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairCasimir::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  for (int i = 1; i < np1; i++)
    for (int j = i; j < np1; j++) setflag[i][j] = 0;

  memory->create(cutsq, np1, np1, "pair:cutsq");
  memory->create(cut, np1, np1, "pair:cut");
  memory->create(c6, np1, np1, "pair:c6");
  memory->create(lam, np1, np1, "pair:lam");
  memory->create(offset, np1, np1, "pair:offset");
}

/* ----------------------------------------------------------------------
   pair_style casimir <global_cutoff>
------------------------------------------------------------------------- */

void PairCasimir::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Illegal pair_style casimir command");
  cut_global = utils::numeric(FLERR, arg[0], false, lmp);

  if (allocated) {
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   pair_coeff I J <C6> <lambda> [cutoff]
------------------------------------------------------------------------- */

void PairCasimir::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR, "Incorrect args for pair_coeff casimir");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double c6_one = utils::numeric(FLERR, arg[2], false, lmp);
  double lam_one = utils::numeric(FLERR, arg[3], false, lmp);
  if (lam_one <= 0.0) error->all(FLERR, "casimir lambda must be positive");

  double cut_one = cut_global;
  if (narg == 5) cut_one = utils::numeric(FLERR, arg[4], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      c6[i][j] = c6_one;
      lam[i][j] = lam_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }
  if (count == 0) error->all(FLERR, "Incorrect args for pair_coeff casimir");
}

/* ---------------------------------------------------------------------- */

double PairCasimir::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    c6[i][j] = 0.0;
    lam[i][j] = 1.0;
    cut[i][j] = cut_global;
  }

  if (offset_flag && cut[i][j] > 0.0) {
    double rc = cut[i][j];
    offset[i][j] = -c6[i][j] / pow(rc, 6.0) / (1.0 + rc / lam[i][j]);
  } else {
    offset[i][j] = 0.0;
  }

  c6[j][i] = c6[i][j];
  lam[j][i] = lam[i][j];
  offset[j][i] = offset[i][j];
  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

void PairCasimir::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, r, rinv, r6inv, evdwl, fpair;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq >= cutsq[itype][jtype]) continue;

      r = sqrt(rsq);
      rinv = 1.0 / r;
      r6inv = rinv * rinv * rinv * rinv * rinv * rinv;

      /* ================= YOUR PHYSICS GOES HERE ======================
         E(r)     = -C6 * r^-6 * S(r)
         S(r)     =  1 / (1 + r/lambda)
         dS/dr    = -1 / (lambda * (1 + r/lambda)^2)
         dE/dr    =  C6 * ( 6 r^-7 S  -  r^-6 dS/dr )

         LAMMPS convention: fpair = -(1/r) dE/dr, and the force on atom i
         is fpair * (x_i - x_j).                                        */

      double u = 1.0 + r / lam[itype][jtype];
      double s = 1.0 / u;
      double ds = -1.0 / (lam[itype][jtype] * u * u);

      double e = -c6[itype][jtype] * r6inv * s;
      double dedr = c6[itype][jtype] * (6.0 * r6inv * rinv * s - r6inv * ds);

      fpair = -dedr * rinv;
      evdwl = e - offset[itype][jtype];

      /* ============== END OF USER-MODIFIABLE BLOCK =================== */

      f[i][0] += delx * fpair;
      f[i][1] += dely * fpair;
      f[i][2] += delz * fpair;
      if (newton_pair || j < nlocal) {
        f[j][0] -= delx * fpair;
        f[j][1] -= dely * fpair;
        f[j][2] -= delz * fpair;
      }

      if (evflag)
        ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   needed by compute group/group, pair_write, and fix numdiff
------------------------------------------------------------------------- */

double PairCasimir::single(int /*i*/, int /*j*/, int itype, int jtype,
                           double rsq, double /*factor_coul*/,
                           double factor_lj, double &fforce)
{
  double r = sqrt(rsq);
  double rinv = 1.0 / r;
  double r6inv = rinv * rinv * rinv * rinv * rinv * rinv;

  double u = 1.0 + r / lam[itype][jtype];
  double s = 1.0 / u;
  double ds = -1.0 / (lam[itype][jtype] * u * u);

  double e = -c6[itype][jtype] * r6inv * s;
  double dedr = c6[itype][jtype] * (6.0 * r6inv * rinv * s - r6inv * ds);

  fforce = factor_lj * (-dedr * rinv);
  return factor_lj * (e - offset[itype][jtype]);
}
