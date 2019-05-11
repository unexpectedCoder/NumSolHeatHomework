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
    double t0 = 80.0;       // Cel degrees
    double ta = 20.0;       // ...........
    double delta_t = 5.0;
    double epsilon = 0.95;  // Blackness
    double c = 40.0;        // Specific heat
    double H = 0.2;
    // For material
    double T[] = { 0.0 + T_ABS, 100.0 + T_ABS, 200.0 + T_ABS };
    double lam[] = { 106.0, 109.0, 110.0 };
    double crho[] = { 3.25e6, 3.34e6, 3.42e6 };

    Wall wall(0.0, 0.1, 10, "steel");
//    wall.setLambda("../NumSolHeatHomework/lambda(T).txt");
    wall.setLambda(T, lam, n1);
    wall.set_crho(T, crho, n1);
    wall.setDens(7800);
    wall.setBlackness(epsilon);
    wall.setSpecificHeat(c);

    Walls walls;
    walls.push_back(wall);

    StartConds sc(t0);
    sc.setGeometry(walls, H);

    // Solution
    ImplicitDiffSchemeCyl solver;
    solver.addWall(wall);
    solver.setStartConds(sc);
    solver.setEnvironment(ta, "../NumSolHeatHomework/env_data.txt");
    solver.solve(dt, delta_t);
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
