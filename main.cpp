#include <QApplication>

#include <iostream>
#include <string>

#include "types.h"


using namespace std;

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);

  try
  {
    Wall wall(0.0, 0.1, 1000);

    wall.setLambda("/home/vadim/QtProjects/NumSolHeatHomework/lambda(T).txt");
//    double T[] = { 0, 100.0, 200.0 };
//    double lam[] = { 1.0, 2.1, 3.5 };
//    wall.setLambda(T, lam, 3);
  }
  catch (const string& ex)
  {
    cout << ex;
  }
  cout << "\nOK!\n\n";

  return a.exec();
}
