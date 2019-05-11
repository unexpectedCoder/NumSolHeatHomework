#include "implicit_diff_scheme_cyl.h"

#include <iostream>
#include <math.h>


using namespace std;


ImplicitDiffSchemeCyl::ImplicitDiffSchemeCyl() :
  is_walls(false), is_startConds(false),
  is_bound1(false), is_bound2(false), is_env(false),
  curTotalN(0), wallsN(0), t_ind(0)
{
  lInterp = gsl_interp_eval;
  sInterp = gsl_spline_eval;
}


ImplicitDiffSchemeCyl::~ImplicitDiffSchemeCyl()
{
  for (size_t i = 0; i < wallsN; ++i)
  {
    delete [] a[i];
    delete [] b[i];
  }
  delete [] a;
  delete [] b;

  // Clear interpolation
//  for (size_t i = 0; i < wallsN; ++i)
//    gsl_spline_free(lLam[i]);
//  delete [] lLam;

  gsl_spline_free(sEnv_c);
  gsl_spline_free(sEnv_Pr);
  gsl_spline_free(sEnv_mu);
  gsl_spline_free(sEnv_nu);
  gsl_spline_free(sEnv_lam);
  gsl_spline_free(sEnv_rho);
  gsl_spline_free(sEnv_a);

  gsl_interp_accel_free(acc);
}


void ImplicitDiffSchemeCyl::addWall(const Wall &w)
{
  walls.push_back(w);
  wallsN++;
  is_walls = true;
}


void ImplicitDiffSchemeCyl::setStartConds(const StartConds &sc)
{
  time = sc.time;
  H = sc.H;
  D = sc.D;

  T0 = sc.T0;
  Tl = T0;
  Tc = T0;
  Tr = T0;

  Tw_vec.push_back(T0);
  t_vec.push_back(time);

  setStartTemperature();

  is_startConds = true;
}


void ImplicitDiffSchemeCyl::setFirstBound(const BoundCond &bc)
{
  if (bc.type != 2)
    throw err.sendEx("this condition type is not supported for the first bound");
  bound1 = bc;
  is_bound1 = true;
}


void ImplicitDiffSchemeCyl::setSecondBound(const BoundCond &bc)
{
  if (bc.type != 3)
    throw err.sendEx("this condition type is not supported for the first bound");
  bound2 = bc;
  is_bound2 = true;
}


void ImplicitDiffSchemeCyl::setEnvironment(double t_amb_C, const string &src_path)
{
  if (t_amb_C < -T_ABS)
    throw err.sendEx("temperature is set less than absolute 0");

  fstream file(src_path.c_str(), ios_base::in);
  if (!file.is_open())
    throw err.sendEx("file is not opened");
  file.close();

  env.Ta = t_amb_C + T_ABS;
  env.dataSize = calcEnvSize(src_path);
  giveMemEnv();
  readEnvData(src_path);

  is_env = true;

  cout << env << '\n';
}


void ImplicitDiffSchemeCyl::solve(double dt, double delta_T)
{
  // Flags checking
  if (!is_env)
    throw err.sendEx("environment is not initialized");
  if (!is_walls)
    throw err.sendEx("walls are not initialized");
  if (!is_bound1)
    throw err.sendEx("first boundary condition is not initialized");
  if (!is_bound2)
    throw err.sendEx("second boundary condition is not initialized");
  if (!is_startConds)
    throw err.sendEx("start conditions are not initialized");
  // Input checking
  if (delta_T < 0.0)
    throw err.sendEx("temperature is set less than ambient temperature");
  if (dt < 0.0)
    throw err.sendEx("time step must be > 0");

  double T_end = env.Ta + delta_T;

  acc = gsl_interp_accel_alloc();
  giveMemDF();
  prepareLInterp();
  prepareSInterp();
  calcDF();
}


void ImplicitDiffSchemeCyl::showWalls() const
{
  if (!walls.empty())
  {
    cout << "Amount of walls: " << wallsN << '\n';
    for (WallCItr i = walls.begin(); i != walls.end(); ++i)
      cout << *i;
  }
}


// *** PRIVATE ***
void ImplicitDiffSchemeCyl::setStartTemperature()
{
  if (T0 < 0.0)
    throw err.sendEx("invalid temperature (less than absolute 0)");

  for (WallItr itr = walls.begin(); itr != walls.end(); ++itr)
  {
    if (itr->is_T)
      throw err.sendEx("T(tau) was already initialized");

    for (size_t i = 0; i < itr->N; ++i)
    {
      itr->r[0] = itr->step * i;
      itr->T[i] = T0;
    }
    itr->is_T = true;
  }
}


size_t ImplicitDiffSchemeCyl::calcEnvSize(const string &path)
{
  size_t s = 0;
  fstream f(path.c_str(), ios_base::in);
  string buf;

  while (!f.eof())
  {
    for (int i = 0; i < 8; ++i)
      f >> buf;
    s++;
  }

  f.close();
  return --s;
}


void ImplicitDiffSchemeCyl::giveMemEnv()
{
  env.T = new double[env.dataSize];
  env.a = new double[env.dataSize];
  env.c = new double[env.dataSize];
  env.Pr = new double[env.dataSize];
  env.mu = new double[env.dataSize];
  env.nu = new double[env.dataSize];
  env.rho = new double[env.dataSize];
  env.lambda = new double[env.dataSize];
}


void ImplicitDiffSchemeCyl::readEnvData(const string &path)
{
  fstream f(path.c_str(), ios_base::in);

  for (size_t i = 0; i < env.dataSize; ++i)
  {
    f >> env.T[i] >> env.lambda[i]
      >> env.rho[i] >> env.c[i]
      >> env.a[i] >> env.nu[i]
      >> env.mu[i] >> env.Pr[i];
  }

  f.close();
}


void ImplicitDiffSchemeCyl::prepareLInterp()
{
  // Function uses linear interpolation for materials
  // because temperature changes in narrow interval

  lLam = new gsl_spline*[wallsN];
  l_c = new gsl_spline*[wallsN];
  for (size_t i = 0; i < wallsN; ++i)
  {
    // lambda
    lLam[i] = gsl_spline_alloc(gsl_interp_cspline, walls[i].dataSize);
    gsl_spline_init(lLam[i],
                    walls[i].T_table, walls[i].lambda,
                    walls[i].dataSize);
    // c
    l_c[i] = gsl_spline_alloc(gsl_interp_cspline, walls[i].dataSize);
    gsl_spline_init(l_c[i],
                    walls[i].T_table, walls[i].c,
                    walls[i].dataSize);
  }
  double a = 2;
}


void ImplicitDiffSchemeCyl::prepareSInterp()
{
  sEnv_c = gsl_spline_alloc(gsl_interp_cspline, env.dataSize);
  sEnv_Pr = gsl_spline_alloc(gsl_interp_cspline, env.dataSize);
  sEnv_mu = gsl_spline_alloc(gsl_interp_cspline, env.dataSize);
  sEnv_nu = gsl_spline_alloc(gsl_interp_cspline, env.dataSize);
  sEnv_lam = gsl_spline_alloc(gsl_interp_cspline, env.dataSize);
  sEnv_rho = gsl_spline_alloc(gsl_interp_cspline, env.dataSize);
  sEnv_a = gsl_spline_alloc(gsl_interp_cspline, env.dataSize);

  gsl_spline_init(sEnv_c, env.T, env.c, env.dataSize);
  gsl_spline_init(sEnv_Pr, env.T, env.Pr, env.dataSize);
  gsl_spline_init(sEnv_mu, env.T, env.mu, env.dataSize);
  gsl_spline_init(sEnv_nu, env.T, env.nu, env.dataSize);
  gsl_spline_init(sEnv_lam, env.T, env.lambda, env.dataSize);
  gsl_spline_init(sEnv_rho, env.T, env.rho, env.dataSize);
  gsl_spline_init(sEnv_a, env.T, env.a, env.dataSize);
}


void ImplicitDiffSchemeCyl::giveMemDF()
{
  a = new double*[wallsN];
  A = new double*[wallsN];
  b = new double*[wallsN];
  B = new double*[wallsN];
  for (size_t i = 0; i < wallsN; ++i)
  {
    a[i] = new double[walls[i].N];
    A[i] = new double[walls[i].N];
    b[i] = new double[walls[i].N];
    B[i] = new double[walls[i].N];
  }
}


void ImplicitDiffSchemeCyl::calcDF()
{
  for (size_t wi = 0; wi < wallsN; ++wi)
  {
    setStartDF(wi);
    for (size_t i = 1; i < walls[wi].N - 1; ++i)
    {
      double *buf = calcTempCoeff(wi, i, false);
      double a1 = buf[0];
      double a2 = buf[1];
    }
  }
}


void ImplicitDiffSchemeCyl::setStartDF(size_t i)
{
  if (i > 0)
  {
    a[i][0] = a[i - 1][walls[i - 1].N - 1];
    b[i][0] = b[i - 1][walls[i - 1].N - 1];
    return;
  }
  if (fabs(bound1.q - 0.0) < EPS)           // i == 0
  {
    a[0][0] = 1.0;
    b[0][0] = 0.0;
    return;
  }
  throw err.sendEx("program is not ready for this shit yet");
}


double* ImplicitDiffSchemeCyl::calcTempCoeff(size_t wi, size_t i, bool is_joint)
{
  if (is_joint && wi == wallsN - 1)
    throw err.sendEx("the last wall doesn't have outer joint point");

  double *res = new double[2];
//  double theta_p = 0.5 * (walls[wi].T[i] + walls[wi].T[i + 1]);
//  double theta_m = 0.5 * (walls[wi].T[i - 1] + walls[wi].T[i]);

//  if (is_joint)
//  {
//    double c1 = lInterp(l_c[wi], walls[wi].T_table, walls[wi].c, walls[wi].T[i], acc);
//    double c2 = lInterp(l_c[wi + 1], walls[wi + 1].T_table, walls[wi + 1].c, walls[wi + 1].T[0], acc);
//    double rho1 = lInterp(l_rho[wi], walls[wi].T_table, walls[wi].rho, , acc);
//    double rho2 = lInterp(l_rho[wi + 1], walls[wi + 1].T_table, walls[wi + 1].rho, , acc);
//    double crho_ = (c1 * rho1 * walls[wi].step + c2 * rho2 * walls[wi + 1].step)
//                   / (walls[wi].step + walls[wi + 1].step);

//    double lam1 = lInterp(lLam[wi], walls[wi].T_table)
//  }

  double theta_p = 0.5 * (walls[wi].T[i] + walls[wi].T[i + 1]);
  double theta_m = 0.5 * (walls[wi].T[i - 1] + walls[wi].T[i]);

  double c = sInterp(l_c[wi], walls[wi].T[i], acc);
  double rho = walls[wi].rho;

  double lam1 = sInterp(lLam[wi], theta_p, acc);
  double lam2 = sInterp(lLam[wi], theta_m, acc);

  res[0] = lam1 / (c * rho);
  res[1] = lam2 / (c * rho);

  return res;
}
