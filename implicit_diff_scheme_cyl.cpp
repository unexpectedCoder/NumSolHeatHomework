#include "implicit_diff_scheme_cyl.h"

#include <iostream>
#include <math.h>


using namespace std;


ImplicitDiffSchemeCyl::ImplicitDiffSchemeCyl() :
  is_walls(false), is_startConds(false),
  is_bound1(false), is_bound2(false), is_env(false),
  totalN(0), wallsN(0), t_ind(0), alphaS(0.0)
{
  lInterp = gsl_interp_eval;
  sInterp = gsl_spline_eval;
}


ImplicitDiffSchemeCyl::~ImplicitDiffSchemeCyl()
{
  delete [] a;
  delete [] b;

  // Clear interpolation
  for (size_t i = 0; i < wallsN; ++i)
  {
    gsl_spline_free(lLam[i]);
    gsl_spline_free(l_c[i]);
  }
  delete [] lLam;
  delete [] l_c;

  gsl_spline_free(sEnv_c);
  gsl_spline_free(sEnv_Pr);
  gsl_spline_free(sEnv_mu);
  gsl_spline_free(sEnv_nu);
  gsl_spline_free(sEnv_lam);
  gsl_spline_free(sEnv_rho);
  gsl_spline_free(sEnv_a);

  gsl_interp_accel_free(acc);
}


void ImplicitDiffSchemeCyl::setWalls(const Walls &ws)
{
  for (WallCItr i = ws.begin(); i != ws.end(); ++i)
  {
    walls.push_back(*i);
    wallsN++;
    totalN += i->N;
  }
  totalN--;
  is_walls = true;
}


void ImplicitDiffSchemeCyl::setStartConds(const StartConds &sc)
{
  time = sc.time;
  H = sc.H;
  D = sc.D;

  T0 = sc.T0;

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

  cout << "N_total = " << totalN << '\n';

  acc = gsl_interp_accel_alloc();
  giveMemDF();
  prepareLInterp();
  prepareSInterp();
  calcDF(dt);

  // TODO: make offsets for DF
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
  // Start temperature destribution

  if (T0 < 0.0)
    throw err.sendEx("invalid temperature (less than absolute 0)");

  theta_buf.resize(totalN);
  for (size_t i = 0; i < totalN; ++i)
    theta_buf[i] = T0;
  theta.push_back(theta_buf);

  for (WallItr itr = walls.begin(); itr != walls.end(); ++itr)
  {
    if (itr->is_T)
      throw err.sendEx("T(tau) was already initialized");

    for (size_t i = 0; i < itr->N; ++i)
      itr->T[i] = T0;
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
  a = new double[totalN];
  A = new double[totalN];
  b = new double[totalN];
  B = new double[totalN];
}


void ImplicitDiffSchemeCyl::calcDF(double dt)
{
  setStartDF();

  size_t offset = 0;
  size_t wi = 0;
  size_t n = walls[wi].N - 1;

  for (size_t i = 1; i < totalN - 1; ++i)
  {
    if (i == n)
    {
      // Joint place

      cout << "i = " << i << '\n';
      cout << "n = " << n << '\n';
      cout << "offset = " << offset << "\n\n";

      if (wi != wallsN - 1)
      {
//        TODO: calculate joint DF
      }

      offset += walls[wi].N - wi;
      n += walls[wi + 1].N - 1;
      wi++;

      continue;
    }

//    double *buf = calcTempCoeff(wi, i - (offset - 1), false);
//    double a1 = buf[0];
//    double a2 = buf[1];

//    A[i] = a1 * dt * (1.0 + 1.0 / (2.0 * i)) / pow(walls[wi].step, 2.0);
//    B[i] = a2 * dt * (1.0 - 1.0 / (2.0 * i)) / pow(walls[wi].step, 2.0);

//    a[i] = A[i] / (1.0 + A[i] + B[i] * (1.0 - a[i - 1]));
//    b[i] = theta[t_ind][i] / A[i]
//           + B[i] / A[i] * a[i - 1] * b[i - 1];
  }
}


void ImplicitDiffSchemeCyl::setStartDF()
{
//  if (i > 0)
//  {
//    a[i][0] = a[i - 1][walls[i - 1].N - 1];
//    b[i][0] = b[i - 1][walls[i - 1].N - 1];
//    A[i][0] = A[i - 1][walls[i - 1].N - 1];
//    B[i][0] = B[i - 1][walls[i - 1].N - 1];
//    return;
//  }

  // i == 0
  if (fabs(bound1.q - 0.0) < EPS)
  {
    a[0] = 1.0;
    b[0] = 0.0;
    A[0] = 0.0;
    B[0] = 0.0;
    return;
  }
  throw err.sendEx("program is not ready for this shit yet");
}


double* ImplicitDiffSchemeCyl::calcTempCoeff(size_t wi, size_t i, bool is_joint)
{
  if (is_joint && wi == wallsN - 1)
    throw err.sendEx("the last wall doesn't have outer joint point");

  double *res = new double[2];
  double theta_p = 0.5 * (walls[wi].T[i] + walls[wi].T[i + 1]);
  double theta_m = 0.5 * (walls[wi].T[i - 1] + walls[wi].T[i]);

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

  double c = sInterp(l_c[wi], walls[wi].T[i], acc);
  double rho = walls[wi].rho;

  cout << theta_p << '\n' /*<< theta_m << '\n'*/;
  double lam1 = sInterp(lLam[wi], theta_p, acc);
  double lam2 = sInterp(lLam[wi], theta_m, acc);

  res[0] = lam1 / (c * rho);
  res[1] = lam2 / (c * rho);

  return res;
}


void ImplicitDiffSchemeCyl::calcTemperature()
{
  calcAlphaSum(theta[t_ind][totalN - 1]);
  for (size_t wi = wallsN - 1; wi + 1 >= 1; --wi)
  {
    double c1 = env.Ta * alphaS * walls[wi].step
                / (sInterp(lLam[wi], theta[t_ind][totalN - 1], acc));
//    double c2 = a[wi][totalN]
  }
}


void ImplicitDiffSchemeCyl::calcAlphaSum(double th)
{
  double T = 0.5 * (th + env.Ta);
  double Gr = g * (th - env.Ta) / (T + env.Ta) * pow(H, 3.0)
              / pow(sInterp(sEnv_nu, T, acc), 2.0);
  double Nu = 0.55 * pow(Gr * sInterp(sEnv_Pr, T, acc), 0.25);
  double al_c = sInterp(sEnv_lam, T, acc) * Nu / H;
  double q_r = C * walls[wallsN - 1].epsilon
               * (pow(th / 100.0, 4.0) - pow(env.Ta / 100.0, 4.0));
  double al_r = q_r / (th - env.Ta);

  alphaS = al_c + al_r;
}
