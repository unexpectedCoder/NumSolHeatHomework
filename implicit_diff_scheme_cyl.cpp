#include "implicit_diff_scheme_cyl.h"

#include <iostream>


using namespace std;


ImplicitDiffSchemeCyl::ImplicitDiffSchemeCyl() :
  wallsNum(0)
{
}


ImplicitDiffSchemeCyl::~ImplicitDiffSchemeCyl()
{
}


void ImplicitDiffSchemeCyl::addWall(const Wall &w)
{
  walls.push_back(w);
  wallsNum++;
}


void ImplicitDiffSchemeCyl::showWalls() const
{
  if (!walls.empty())
  {
    cout << "Amount of walls: " << wallsNum << '\n';
    for (WallCItr i = walls.begin(); i != walls.end(); ++i)
      cout << *i;
  }
}
