#ifndef IMPLICIT_DIFF_SCHEME_CYL_H
#define IMPLICIT_DIFF_SCHEME_CYL_H

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "types.h"

#define RES_PATH "../NumSolHeatHomework/results.txt"


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
  gsl_interp_accel *lAcc;
  gsl_spline **lLam;          // 'l*' means 'linear'
  gsl_spline **l_c;

  double (*lInterp)(const gsl_spline*, double, gsl_interp_accel*);
  // **-pointers are used for each wall

  gsl_interp_accel *sAcc;
  gsl_spline *sEnv_lam;       // 's*' means 'spline'
  gsl_spline *sEnv_rho;
  gsl_spline *sEnv_c;
  gsl_spline *sEnv_a;
  gsl_spline *sEnv_nu;
  gsl_spline *sEnv_mu;
  gsl_spline *sEnv_Pr;

  // Pointer to spline interpolation function
  double (*sInterp)(const gsl_spline*, double, gsl_interp_accel*);
  // .............................................................

  // Others
  size_t t_ind;   // Current time layer index
  double alphaS;  // Summary heat emission coeff
  double *r;      // Common coordinates

  // For the results
  std::vector<double> time_vec;             // Time vector
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
  void setCommonCoords();
  size_t calcEnvSize(const std::string &path);
  void readEnvData(const std::string &path);
  void giveMemEnv();
  void prepareLInterp();
  void prepareSInterp();
  void giveMemDF();
  void calcDF(double dt);
  void calcJointDF(double dt, size_t wi, size_t i);
  void calcInnerDF(double dt, size_t wi, size_t i);
  void setStartDF();
  double* calcTempCoeffs(size_t wi, size_t i);
  double calcJointTempCoeff(size_t wi, size_t i);
  void calcTemperature();
  void calcAlphaSum(double th);
  void writeResultsFile(const std::string &path);
};


#endif // IMPLICIT_DIFF_SCHEME_CYL_H
