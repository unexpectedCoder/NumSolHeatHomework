#ifndef IMPLICIT_DIFF_SCHEME_CYL_H
#define IMPLICIT_DIFF_SCHEME_CYL_H

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "types.h"


class ImplicitDiffSchemeCyl
{
private:
  Error err;

  // Flags that solver is ready to solve
  bool  is_walls,
        is_startConds,
        is_bound1, is_bound2,
        is_env;

  size_t totalN;

  // For init
  Walls walls;                // Vector of walls
  size_t wallsN;              // Amount of walls
  BoundCond bound1, bound2;   // Left & right boundary conditions
  Environment env;            // Environment data
  double H, D;                // Outer geometry
  double T0;                  // Start temperature
  double time;                // Current time

  // Driving factors (DF)
  double *a, *A;
  double *b, *B;

  // For interpolation
  gsl_interp_accel *acc;
  gsl_spline **lLam;          // 'l*' means 'linear'
  gsl_spline **l_c;
  // **-pointers are used for each wall
  double (*lInterp)(const gsl_interp*, const double*, const double*, double, gsl_interp_accel*);

  gsl_spline *sEnv_lam;       // 's*' means 'spline'
  gsl_spline *sEnv_rho;
  gsl_spline *sEnv_c;
  gsl_spline *sEnv_a;
  gsl_spline *sEnv_nu;
  gsl_spline *sEnv_mu;
  gsl_spline *sEnv_Pr;
  double (*sInterp)(const gsl_spline*, double, gsl_interp_accel*);
  // .................

  size_t t_ind;               // Current time layer index
  double alphaS;

  // For the results
  std::vector<double> t_vec;                // Time vector
  std::vector<std::vector<double> > theta;  // Wall inner temperature field
  std::vector<double> theta_buf;            // 1D vector for the adding to 2D theta vector
  std::vector<double> Tw_vec;               // Wall outer temperature vector

public:
  ImplicitDiffSchemeCyl();
  ~ImplicitDiffSchemeCyl();

  void setWalls(const Walls &ws);
  void setStartConds(const StartConds &sc);
  void setFirstBound(const BoundCond &bc);
  void setSecondBound(const BoundCond &bc);
  void setEnvironment(double t_amb_C, const std::string &src_path);

  void solve(double dt, double t_end_C);

  // Out funcs
  void showWalls() const;

private:
  void setStartTemperature();
  size_t calcEnvSize(const std::string &path);
  void readEnvData(const std::string &path);
  void giveMemEnv();
  void prepareLInterp();
  void prepareSInterp();
  void giveMemDF();
  void calcDF(double dt);
  void setStartDF();
  double* calcTempCoeff(size_t wi, size_t i, bool is_joint);
  void calcTemperature();
  void calcAlphaSum(double th);
};


#endif // IMPLICIT_DIFF_SCHEME_CYL_H
