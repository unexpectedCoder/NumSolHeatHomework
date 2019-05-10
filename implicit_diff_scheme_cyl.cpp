#include "implicit_diff_scheme_cyl.h"

#include <iostream>
#include <math.h>


using namespace std;


ImplicitDiffSchemeCyl::ImplicitDiffSchemeCyl() :
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
}


void ImplicitDiffSchemeCyl::setStartConds(const StartConds &sc)
{
  H = sc.H;
  D = sc.D;
  T0 = sc.T;
  T_amb = sc.T_amb;
  time = sc.time;
}


void ImplicitDiffSchemeCyl::setFirstBound(const BoundConds &bc)
{
  if (bc.type != 2)
    throw err.sendEx("this condition type is not supported for the first bound");
  bound1 = bc;
}


void ImplicitDiffSchemeCyl::setSecondBound(const BoundConds &bc)
{
  if (bc.type != 3)
    throw err.sendEx("this condition type is not supported for the first bound");
  bound2 = bc;
}


void ImplicitDiffSchemeCyl::solve(double dt, double t_end_C)
{
  if (t_end_C < -T_ABS)
    throw err.sendEx("temperature is set less than absolute 0");
  if (dt < 0.0)
    throw err.sendEx("time step must be > 0");

  double T_end = t_end_C + T_ABS;

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
  b = new double*[wallsN];
  for (size_t i = 0; i < wallsN; ++i)
  {
    a[i] = new double[walls[i].N];
    b[i] = new double[walls[i].N];
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
  if (fabs(bound1.q - 0.0) < EPS)
  {
    a[i][0] = 1.0;
    b[i][0] = 0.0;
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
