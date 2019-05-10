#ifndef IMPLICIT_DIFF_SCHEME_CYL_H
#define IMPLICIT_DIFF_SCHEME_CYL_H

#include <vector>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "types.h"


class ImplicitDiffSchemeCyl
{
private:
  Error err;

  // For init
  Walls walls;
  size_t wallsN;
  BoundConds bound1, bound2;
  double H, D, T0;
  double time;
  double T_amb;

  // Driving factors (DF)
  std::vector<double> a;
  std::vector<double> b;

  // For interpolation
  gsl_interp_accel *acc;
  gsl_interp **lLam;      // 'l*' means 'linear'
  // .................

public:
  ImplicitDiffSchemeCyl();
  ~ImplicitDiffSchemeCyl();

  void addWall(const Wall &w);
  void setStartConds(const StartConds &sc);
  void setFirstBound(const BoundConds &bc);
  void setSecondBound(const BoundConds &bc);

  void solve(double dt, double t_end_C);

  // Out funcs
  void showWalls() const;

private:
  void prepareInterp();
  void calcDF();
  void setStartDF();
  void freeInterp();
};


#endif // IMPLICIT_DIFF_SCHEME_CYL_H
