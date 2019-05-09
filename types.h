#ifndef TYPES_H
#define TYPES_H

#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#define EPS 1e-12


struct Error
{
  const std::string& sendEx(const std::string& ex)
  {
    return ex;
  }
};


// *** WALL ***
struct Wall
{
  Error err;

  int N;
  double r1, r2;
  double step;
  double **lambda_T;
  size_t lamSize;
  bool is_lambda;

  Wall(double r1, double r2, int n) : lamSize(0), is_lambda(false)
  {
    if (r1 < 0.0 || fabs(r2 - r1) < EPS)
      throw err.sendEx("\n\tError: r1 < 0 or r2 < r1!\n");
    if (n < 2)
      throw err.sendEx("\n\tError: number of segments < 1!\n");

    this->r1 = r1;
    this->r2 = r2;
    N = n;
    step = (r2 - r1) / N;

    lambda_T = new double*[2];
  }

  ~Wall()
  {
    for (size_t i = 0; i < 2; ++i)
      delete [] lambda_T[i];
    delete [] lambda_T;
  }

  void setLambda(const std::string& file_path);
  void setLambda(const double *T, const double *lam, size_t n);
};


// ___ WALL FUNCS ___
void Wall::setLambda(const std::string& file_path)
{
  if (is_lambda)
    throw err.sendEx("\n\tError: lamba is already set!\n");
  is_lambda = true;

  std::ifstream file(file_path.c_str(), std::ios_base::out);
  if (!file.is_open())
    throw err.sendEx("\n\tError: file is not opened!\n");

  double buf;
  while (!file.eof())
  {
    file >> buf >> buf;
    ++lamSize;
  }
  lamSize--;
  file.close();

  file.open(file_path.c_str(), std::ios_base::out);
  for (size_t i = 0; i < 2; ++i)
    lambda_T[i] = new double[lamSize];
  for (size_t i = 0; i < lamSize; ++i)
    file >> lambda_T[0][i] >> lambda_T[1][i];
  file.close();
}

void Wall::setLambda(const double *T, const double *lam, size_t n)
{
  if (is_lambda)
    throw err.sendEx("\n\tError: lamba is already set!\n");
  is_lambda = true;
  if (n < 2)
    throw err.sendEx("\n\tError: arrays' size < 2!\n");

  lamSize = n;
  for (size_t i = 0; i < 2; ++i)
    lambda_T[i] = new double[lamSize];

  for (size_t i = 0; i < lamSize; ++i)
  {
    lambda_T[0][i] = T[i];
    lambda_T[1][i] = lam[i];
  }
}


#endif // TYPES_H
