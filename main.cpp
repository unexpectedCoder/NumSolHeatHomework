//#include <QApplication>

#include <iostream>
#include <string>

#include "implicit_diff_scheme_cyl.h"


using namespace std;

int main(int argc, char *argv[])
{
  try
  {
    const size_t n1 = 3;

    // Initialization
    double dt = 20.0;       // seconds
    double t0 = 174.0;      // Cel degrees
    double ta = 20.0;       // ...........
    double delta_t = 5.0;
    double epsilon = 0.95;  // Blackness
    double H = 0.2;
    // For material
    double T_table[] = { 0.0, 100.0, 200.0 };
    double lam[] = { 106.0, 109.0, 110.0 };
    double rho = 8500.0;
    double c[] = { 3.25e6, 3.34e6, 3.42e6 };
    for (size_t i = 0; i < n1; ++i)
      c[i] /= rho;

    Wall wall(0.0, 0.1, 10, "steel");
//    wall.setLambda("../NumSolHeatHomework/lambda(T).txt");
    wall.setLambdaT(T_table, lam, n1);
    wall.setDens(rho);
    wall.setBlackness(epsilon);
    wall.setSpecificHeat(c, n1);

    Walls walls;
    walls.push_back(wall);

    StartConds sc(t0);
    sc.setGeometry(walls, H);

    BoundCond bc1, bc2;
    bc1.setType2(0.0);
    bc2.setType3(ta);

    // Solution
    ImplicitDiffSchemeCyl solver;
    // ...set walls
    solver.addWall(wall);
    // ...for init control
    solver.showWalls();
    // ...set bounds
    solver.setFirstBound(bc1);
    solver.setSecondBound(bc2);
    // ...set start conditions and environment
    solver.setStartConds(sc);
    solver.setEnvironment(ta, "../NumSolHeatHomework/env_data.txt");

    solver.solve(dt, delta_t);

    // ...for init control
    solver.showWalls();
  }
  catch (const string& ex)
  {
    cout << ex;
  }

  return 0;
//  QApplication a(argc, argv);
//  return a.exec();
}
