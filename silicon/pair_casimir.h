/* ----------------------------------------------------------------------
   pair_casimir.h

   Minimal but complete custom pair style for LAMMPS.

   E(r) = -C6 / r^6 * S(r),   S(r) = 1 / (1 + r/lambda)

   S(r) is a crude retardation switch: it turns the non-retarded r^-6
   van der Waals form into the retarded r^-7 Casimir-Polder form for
   r >> lambda.  Replace the marked block in compute() with whatever
   functional form you actually want -- everything else is boilerplate.

   Usage in an input script:
     pair_style  hybrid/overlay sw sw casimir 30.0
     pair_coeff  * * sw 1 Si.sw Si NULL
     pair_coeff  * * sw 2 Si.sw NULL Si
     pair_coeff  1 2 casimir <C6> <lambda>
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(casimir,PairCasimir);
// clang-format on
#else

#ifndef LMP_PAIR_CASIMIR_H
#define LMP_PAIR_CASIMIR_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCasimir : public Pair {
 public:
  PairCasimir(class LAMMPS *);
  ~PairCasimir() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  double single(int, int, int, int, double, double, double, double &) override;

 protected:
  double cut_global;
  double **cut;
  double **c6;       // dispersion coefficient   [energy * length^6]
  double **lam;      // retardation length scale [length]
  double **offset;   // energy shift at the cutoff

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
