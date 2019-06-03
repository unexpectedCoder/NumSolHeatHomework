#include <QApplication>

#include <iostream>
#include <string>

#include "implicit_diff_scheme_cyl.h"
#include "plotter.h"
#include "mainwindow.h"


using namespace std;

int main(int argc, char *argv[])
{
  try
  {
    const size_t n = 3;

    /* ********************************************* */
    /* ***** For comparing with the experiment ***** */
//    size_t N = 200;         // Number of grid segments
//    double dt = 20.0;       // seconds
//    double t0 = 174.0;      // Cel degrees
//    double ta = 22.2;       // ...........
//    double t_end = 50.0;    // Termination condition
//    double delta_t = t_end - ta;
//    double epsilon = 0.95;  // Blackness
//    double H = 0.25;        // Cyl heigh
//    double D1 = 0.0;        // Inner cyl diameter
//    double D2 = 0.01;       // Outer cyl diameter
//    double D3 = 0.02;
//    double D4 = 0.0299;
//    // For material
//    double T_table[] = { 0.0, 100.0, 200.0 };
//    double lam[] = { 106.0, 109.0, 110.0 };
//    double rho = 8500.0;
//    double c[] = { 3.25e6, 3.34e6, 3.42e6 };
//    for (size_t i = 0; i < n; ++i)
//      c[i] /= rho;

//    Wall wall1(D1 / 2.0, D2 / 2.0, N, "steel");
//    wall1.setGrid();
//    wall1.setLambdaT(T_table, lam, n);
//    wall1.setDens(rho);
//    wall1.setBlackness(epsilon);
//    wall1.setSpecificHeat(c, n);

//    Wall wall2(D2 / 2.0, D3 / 2.0, N, "steel");
//    wall2.setGrid();
//    wall2.setLambdaT(T_table, lam, n);
//    wall2.setDens(rho);
//    wall2.setBlackness(epsilon);
//    wall2.setSpecificHeat(c, n);

//    Wall wall3(D3 / 2.0, D4 / 2.0, N, "steel");
//    wall3.setGrid();
//    wall3.setLambdaT(T_table, lam, n);
//    wall3.setDens(rho);
//    wall3.setBlackness(epsilon);
//    wall3.setSpecificHeat(c, n);

//    Walls walls;
//    walls.push_back(wall1);
//    walls.push_back(wall2);
//    walls.push_back(wall3);
    /* ***** END OF Experimental ***** */
    /* ******************************* */

    // Initialization
    size_t N = 300;         // Number of grid segments
    double dt = 20.0;       // seconds
    double t0 = 90.0;       // Cel degrees
    double ta = 20.0;       // ...........
    double t_end = 30.0;    // Termination condition
    double delta_t = t_end - ta;

    // Cylinder geometry
    double H = 0.9;
    double D1 = 0.0;
    double D2 = 0.032;
    double D3 = 0.08;
    double D4 = 0.09;

    double T_table[] = { 0.0, 100.0, 200.0 };
    // Wall #1
    double lam1[] = { 35.0, 36.0, 37.0 };
    double rho1 = 7850.0;
    double c1[] = { 490.0, 496.0, 504.0 };
    double epsilon1 = 0.8;
    // Wall #2
//    double lam2[] = { 0.25, 0.26, 0.27 };
    double lam2[] = { 20.0, 20.1, 20.3 };
    double rho2 = 1080.0;
    double c2[] = { 3.31e3, 3.315e3, 3.2e3 };
    double epsilon2 = 0.9;
    // Wall #3
    double lam3[] = { 35.0, 36.0, 37.0 };
    double rho3 = 7850.0;
    double c3[] = { 490.0, 496.0, 504.0 };
    double epsilon3 = 0.95;

    Wall wall1(D1 / 2.0, D2 / 2.0, N, "20HGSA");
    wall1.setGrid();
    wall1.setLambdaT(T_table, lam1, n);
    wall1.setDens(rho1);
    wall1.setBlackness(epsilon1);
    wall1.setSpecificHeat(c1, n);

    Wall wall2(D2 / 2.0, D3 / 2.0, N, "POJ-70");
    wall2.setGrid();
    wall2.setLambdaT(T_table, lam2, n);
    wall2.setDens(rho2);
    wall2.setBlackness(epsilon2);
    wall2.setSpecificHeat(c2, n);

    Wall wall3(D3 / 2.0, D4 / 2.0, N, "20HGSA");
    wall3.setGrid();
    wall3.setLambdaT(T_table, lam3, n);
    wall3.setDens(rho3);
    wall3.setBlackness(epsilon3);
    wall3.setSpecificHeat(c3, n);

    Walls walls;
    walls.push_back(wall1);
    walls.push_back(wall2);
    walls.push_back(wall3);

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

  // Graphic part
  QApplication a(argc, argv);

  Plotter plot;
  plot.setData("../NumSolHeatHomework/results.txt", true);
  plot.createChart();
  plot.setAxis(0, 30000, 20, 100, 7, 9);
//  plot.setAxis();

  MainWindow w;
  w.setCentralWidget(plot.getChartView());
  w.resize(1200, 800);
  w.show();

  return a.exec();
}
