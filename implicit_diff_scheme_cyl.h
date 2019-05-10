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
  Walls walls;                // Vector of walls
  size_t wallsN;              // Amount of walls
  BoundConds bound1, bound2;  // Left & right boundary conditions
  double H, D;                // Outer geometry
  double T0;                  // Start temperature
  double time;                // Current time
  double T_amb;               // Ambient temperature
  double Tw;                  // Current wall temperature

  // Driving factors (DF)
  double **a;
  double **b;

  // For interpolation
  gsl_interp_accel *acc;
  gsl_interp **lLam;          // 'l*' means 'linear'
  // Interp func for lambda of wall material
  double (*lLam_f)(const gsl_interp*, const double*, const double*, double, gsl_interp_accel*);
  // .................

  size_t t_ind;               // Current time layer index

  // For the results
  std::vector<double> t_vec;  // Time vector
  std::vector<double> Tw_vec; // Wall temperature vector

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
  void giveMemDF();
  void calcDF();
  void setStartDF(size_t i);
  double calcTempCoeff(size_t i, char plus_minus);
  void freeInterp();
};


#endif // IMPLICIT_DIFF_SCHEME_CYL_H
