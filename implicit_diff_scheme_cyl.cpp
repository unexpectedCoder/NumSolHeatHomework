#include "implicit_diff_scheme_cyl.h"

#include <iostream>
#include <math.h>


using namespace std;


ImplicitDiffSchemeCyl::ImplicitDiffSchemeCyl() :
  wallsN(0)
{}


ImplicitDiffSchemeCyl::~ImplicitDiffSchemeCyl()
{}


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
  for (size_t i = 0; i < wallsN; ++i)
  {
    lLam[i] = gsl_interp_alloc(gsl_interp_linear, walls[i].lamSize);
    gsl_interp_init(lLam[i],
                    walls[i].lambda_T[0], walls[i].lambda_T[1],
                    walls[i].lamSize);
  }
}


void ImplicitDiffSchemeCyl::calcDF()
{
  setStartDF();
}


void ImplicitDiffSchemeCyl::setStartDF()
{
  if (fabs(bound1.q - 0.0) < EPS)
  {
    a.push_back(1.0);
    b.push_back(0.0);
    return;
  }
  throw err.sendEx("program is not ready for this shit yet");
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
