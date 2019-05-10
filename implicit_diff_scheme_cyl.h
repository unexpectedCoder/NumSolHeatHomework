#ifndef IMPLICIT_DIFF_SCHEME_CYL_H
#define IMPLICIT_DIFF_SCHEME_CYL_H

#include "types.h"


class ImplicitDiffSchemeCyl
{
private:
  Walls walls;
  size_t wallsNum;

public:
  ImplicitDiffSchemeCyl();
  ~ImplicitDiffSchemeCyl();

  void addWall(const Wall &w);
  void showWalls() const;
};


#endif // IMPLICIT_DIFF_SCHEME_CYL_H
