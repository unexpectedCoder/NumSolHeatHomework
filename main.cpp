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
    size_t N = 300;         // Number of grid segments
    double dt = 20.0;       // seconds
    double t0 = 174.0;      // Cel degrees
    double ta = 22.2;       // ...........
    double t_end = 50.0;    // Termination condition
    double delta_t = t_end - ta;
    double epsilon = 0.95;  // Blackness
    double H = 0.25;        // Cyl heigh
    double D1 = 0.0;        // Inner cyl diameter
    double D2 = 0.0299;     // Outer cyl diameter
    // For material
    double T_table[] = { 0.0, 100.0, 200.0 };
    double lam[] = { 106.0, 109.0, 110.0 };
    double rho = 8500.0;
    double c[] = { 3.25e6, 3.34e6, 3.42e6 };
    for (size_t i = 0; i < n1; ++i)
      c[i] /= rho;

    Wall wall(D1 / 2.0, D2 / 2.0, N, "steel");
//    wall.setLambda("../NumSolHeatHomework/lambda(T).txt");
    wall.setGrid();
    wall.setLambdaT(T_table, lam, n1);
    wall.setDens(rho);
    wall.setBlackness(epsilon);
    wall.setSpecificHeat(c, n1);

    Wall wall2(0.1, 0.2, 20, "steel");
    wall2.setGrid();
    wall2.setLambdaT(T_table, lam, n1);
    wall2.setDens(rho);
    wall2.setBlackness(epsilon);
    wall2.setSpecificHeat(c, n1);

    Wall wall3(0.2, 0.25, 10, "steel");
    wall3.setGrid();
    wall3.setLambdaT(T_table, lam, n1);
    wall3.setDens(rho);
    wall3.setBlackness(epsilon);
    wall3.setSpecificHeat(c, n1);

    Walls walls;
    walls.push_back(wall);
//    walls.push_back(wall2);
//    walls.push_back(wall3);

    StartConds sc(t0);
    sc.setGeometry(walls, H);

    BoundCond bc1, bc2;
    bc1.setType2(0.0);
    bc2.setType3(ta);

    // Solution
    ImplicitDiffSchemeCyl solver;
    // ...set walls
    solver.setWalls(walls);
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
