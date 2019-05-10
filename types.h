#ifndef TYPES_H
#define TYPES_H

#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#define EPS 1e-12       // Accuracy
#define T_ABS 273.15    // Absolute difference between C and K


struct Error
{
  std::string mess = "";
  const std::string& sendEx(const std::string& ex)
  {
    mess = "\n\tError: " + ex + "!\n";
    return mess;
  }
};


// *** WALL ***
struct Wall
{
  Error err;

  int N;              // Number of spacing segments
  double r1, r2;      // Inner and outer radius of cylinder
  double step;        // Space step
  double **lambda_T;  // lambda(T) - wall's heat transfer coeff
  size_t lamSize;     // Size of lambda data
  bool is_lambda;     // Was the lambda(T) initialized
  double **T;         // T(tau): [0][:] - time, [1][:] - temperature
  double epsilon;     // Blackness
  double rho;         // Material density
  double c;           // Specific heat

  Wall(double r1, double r2, int n) :
    lamSize(0), is_lambda(false), epsilon(1.0)
  {
    if (r1 < 0.0 || fabs(r2 - r1) < EPS)
      throw err.sendEx("r1 < 0 or r2 < r1");
    if (n < 2)
      throw err.sendEx("number of segments < 1");

    this->r1 = r1;
    this->r2 = r2;
    N = n;
    step = (r2 - r1) / N;

    lambda_T = new double*[2];
    T = new double*[2];
  }

  ~Wall()
  {
    for (size_t i = 0; i < 2; ++i)
    {
      delete [] lambda_T[i];
      delete [] T[i];
    }
    delete [] lambda_T;
    delete [] T;
  }

  void setLambda(const std::string& file_path);
  void setLambda(const double *T, const double *lam, size_t n);
  void setStartTemperature(double temp_C);
  void setBlackness(double epsilon);
  void setDens(double rho);
  void setSpecificHeat(double c);
};


// *** Boundary conditions (type 2) ***
struct BoundConds2
{
  int type;
  double q;

  BoundConds2(double q) : type(2), q(q) {}
  ~BoundConds2() {}
};


// *** Boundary conditions (type 3) ***
struct BoundConds3
{
  int type;
  double T_amb, alpha;

  BoundConds3(double t_ambient_C, double alpha) :
    type(3), T_amb(t_ambient_C + T_ABS), alpha(alpha) {}
  ~BoundConds3() {}
};


#endif // TYPES_H
