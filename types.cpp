#include "types.h"


using namespace std;


// ___ WALL FUNCS ___
void Wall::setLambda(const string& file_path)
{
  /*
   * Read function lambda(T) from source file
   * and set table **<var>[2][size] with T and lambda.
  */

  if (is_lambda)
    throw err.sendEx("lamba is already set!");
  is_lambda = true;

  ifstream file(file_path.c_str(), ios_base::out);
  if (!file.is_open())
    throw err.sendEx("file is not opened");

  double buf;
  while (!file.eof())
  {
    file >> buf >> buf;
    ++lamSize;
  }
  lamSize--;
  file.close();

  file.open(file_path.c_str(), ios_base::out);
  for (size_t i = 0; i < 2; ++i)
    lambda_T[i] = new double[lamSize];
  for (size_t i = 0; i < lamSize; ++i)
  {
    file >> lambda_T[0][i] >> lambda_T[1][i];
    lambda_T[0][i] += T_ABS;
  }
  file.close();
}


void Wall::setLambda(const double *T, const double *lam, size_t n)
{
  /*
   * Set table **<var>[2][size] with T and lambda
   * using source array T and lambda.
  */

  if (is_lambda)
    throw err.sendEx("lamba is already set");
  is_lambda = true;
  if (n < 2)
    throw err.sendEx("arrays' size < 2");

  lamSize = n;
  for (size_t i = 0; i < 2; ++i)
    lambda_T[i] = new double[lamSize];
  for (size_t i = 0; i < lamSize; ++i)
  {
    lambda_T[0][i] = T[i] + T_ABS;
    lambda_T[1][i] = lam[i];
  }
}


void Wall::setStartTemperature(double temp_C)
{
  if (temp_C + T_ABS < 0.0)
    throw err.sendEx("invalid temperature (less than absolute 0)");

  for (int i = 0; i < 2; ++i)
    T[i] = new double[N];
  for (int i = 0; i < N; ++i)
  {
    T[0][i] = 0.0;
    T[1][i] = temp_C + T_ABS;
  }
}


void Wall::setBlackness(double epsilon)
{
  if (epsilon < 0.0 || epsilon > 1.0)
    throw err.sendEx("blackness must be in [0; 1]");
  this->epsilon = epsilon;
}


void Wall::setDens(double rho)
{
  if (rho < 0.0)
    throw err.sendEx("density must be greater than 0");
  this->rho = rho;
}


void Wall::setSpecificHeat(double c)
{
  if (c < 0.0)
    throw err.sendEx("specific heat must be greater than 0");
  this->c = c;
}
