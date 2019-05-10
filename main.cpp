//#include <QApplication>

#include <iostream>
#include <string>

#include "implicit_diff_scheme_cyl.h"


using namespace std;

int main(int argc, char *argv[])
{
  try
  {
    StartConds sc(80.0, 15.0, 0.1, 0.2);

    Wall wall(0.0, 0.1, 10, "steel");

//    wall.setLambda("../NumSolHeatHomework/lambda(T).txt");
    double T[] = { 0.0 + T_ABS, 100.0 + T_ABS, 200.0 + T_ABS };
    double lam[] = { 106.0, 109.0, 110.0 };
    wall.setLambda(T, lam, 3);
    wall.setStartTemperature(80.0);
    wall.setDens(7800);
    wall.setBlackness(0.5);
    wall.setSpecificHeat(40);
    cout << wall << endl;

    ImplicitDiffSchemeCyl solver;
    solver.addWall(wall);
    solver.solve(20.0, sc.T_amb);
    solver.showWalls();
  }
  catch (const string& ex)
  {
    cout << ex;
  }
  cout << "\nOK!\n\n";

  return 0;
//  QApplication a(argc, argv);
//  return a.exec();
}
