#include "implicit_diff_scheme_cyl.h"

#include <iostream>
#include <math.h>


using namespace std;


ImplicitDiffSchemeCyl::ImplicitDiffSchemeCyl() :
  is_walls(false), is_startConds(false),
  is_bound1(false), is_bound2(false), is_env(false),
  wallsN(0), t_ind(0)
{
  interp_f = gsl_interp_eval;
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
}


void ImplicitDiffSchemeCyl::addWall(const Wall &w)
{
  walls.push_back(w);
  wallsN++;
  is_walls = true;
}


void ImplicitDiffSchemeCyl::setStartConds(const StartConds &sc)
{
  H = sc.H;
  D = sc.D;
  T0 = sc.T0;
  time = sc.time;

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

  prepareInterp();
  giveMemDF();
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
      itr->T[0][i] = itr->step * i;
      itr->T[1][i] = T0;
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


void ImplicitDiffSchemeCyl::prepareInterp()
{
  acc = gsl_interp_accel_alloc();

  // Using linear interpolation for materials
  // because temperature changes in narrow interval
  lLam = new gsl_interp*[wallsN];
  l_crho = new gsl_interp*[wallsN];
  for (size_t i = 0; i < wallsN; ++i)
  {
    // Lambda
    lLam[i] = gsl_interp_alloc(gsl_interp_linear, walls[i].dataSize);
    gsl_interp_init(lLam[i],
                    walls[i].lambda_T[0], walls[i].lambda_T[1],
                    walls[i].dataSize);
    // c*rho
    l_crho[i] = gsl_interp_alloc(gsl_interp_linear, walls[i].dataSize);
    gsl_interp_init(l_crho[i],
                    walls[i].crho[0], walls[i].crho[1],
                    walls[i].dataSize);
  }
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
  for (size_t i = 0; i < wallsN; ++i)
  {
    setStartDF(i);
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


double ImplicitDiffSchemeCyl::calcTempCoeff(size_t i, char plus_minus)
{
  if (plus_minus == '+')
  {
    return -1.0;
  }
  return -1.0;
}


void ImplicitDiffSchemeCyl::freeInterp()
{
  if (lLam)
    for (size_t i = 0; i < wallsN; ++i)
      gsl_interp_free(lLam[i]);
  delete [] lLam;

  if (acc)
    gsl_interp_accel_free(acc);
}
