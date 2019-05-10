#ifndef TYPES_H
#define TYPES_H

#include <fstream>
//#include <iostream>
#include <vector>
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

  Error() {}
  ~Error()
  {
    mess.~basic_string();
  }
};


// *** WALL ***
struct Wall
{
  Error err;

  int N;                  // Number of spacing segments
  double r1, r2;          // Inner and outer radius of cylinder
  double step;            // Space step
  double **lambda_T;      // lambda(T) - wall's heat transfer coeff
  size_t lamSize;         // Size of lambda data
  bool is_lambda;         // Was the lambda(T) initialized
  double **T;             // T(tau): [0][:] - time, [1][:] - temperature
  double epsilon;         // Blackness
  double rho;             // Material density
  double c;               // Specific heat
  std::string material;   // Material's name

  Wall(double r1, double r2, int n, const std::string &material = "NONE") :
    lamSize(0), is_lambda(false), epsilon(1.0), rho(0.0), c(0.0), material(material)
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

  Wall(const Wall &w)
  {
    N = w.N;
    r1 = w.r1;
    r2 = w.r2;
    step = w.step;
    lamSize = w.lamSize;
    is_lambda = w.is_lambda;
    epsilon = w.epsilon;
    rho = w.rho;
    c = w.c;
    material = w.material;

    lambda_T = new double*[2];
    T = new double*[2];
    for (size_t i = 0; i < 2; ++i)
    {
      lambda_T[i] = new double[lamSize];
      T[i] = new double[lamSize];
    }

    for (size_t i = 0; i < lamSize; ++i)
    {
      lambda_T[0][i] = w.lambda_T[0][i];
      lambda_T[1][i] = w.lambda_T[0][i];
      T[0][i] = w.T[0][i];
      T[1][i] = w.T[1][i];
    }
  }

  Wall& operator=(const Wall &w)
  {
    if (this == &w)
      return *this;

    N = w.N;
    r1 = w.r1;
    r2 = w.r2;
    step = w.step;
    lamSize = w.lamSize;
    is_lambda = w.is_lambda;
    epsilon = w.epsilon;
    rho = w.rho;
    c = w.c;
    material = w.material;

    lambda_T = new double*[2];
    T = new double*[2];
    for (size_t i = 0; i < 2; ++i)
    {
      lambda_T[i] = new double[lamSize];
      T[i] = new double[lamSize];
    }

    for (size_t i = 0; i < lamSize; ++i)
    {
      lambda_T[0][i] = w.lambda_T[0][i];
      lambda_T[1][i] = w.lambda_T[0][i];
      T[0][i] = w.T[0][i];
      T[1][i] = w.T[1][i];
    }
    return *this;
  }

  void setLambda(const std::string& file_path);
  void setLambda(const double *T, const double *lam, size_t n);
  void setStartTemperature(double temp_C);
  void setBlackness(double epsilon);
  void setDens(double rho);
  void setSpecificHeat(double c);

  inline friend std::ostream& operator<<(std::ostream &os, const Wall &w);

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
};


std::ostream& operator<<(std::ostream &os, const Wall &w)
{
  os << "\n... Wall ...:\n"
     << "\tmaterial: " << w.material << '\n'
     << "\tN: " << w.N << '\n'
     << "\tr1, m: " << w.r1 << '\n'
     << "\tr2, m: " << w.r2 << '\n'
     << "\tstep (dr), m: " << w.step << '\n'
     << "\trho, kg/m^3: " << w.rho << '\n'
     << "\tc, J/(kg*K): " << w.c << '\n'
     << "\tepsilon: " << w.epsilon << '\n';

  os << "\tlambda(T):\n"
     << "\t\tT, C\tlambda, W/(m*K)\n";
  for (size_t i = 0; i < w.lamSize; ++i)
    os << "\t#" << i + 1 << ".\t"
       << w.lambda_T[0][i] - T_ABS << '\t'
       << w.lambda_T[1][i] << '\n';
  os << '\n';

  return os;
}


typedef std::vector<Wall> Walls;
typedef std::vector<Wall>::iterator WallItr;
typedef std::vector<Wall>::const_iterator WallCItr;
// *** END OF WALL ***


// *** Boundary conditions (type 2) ***
struct BoundConds2
{
  int type;
  double q;

  BoundConds2(double q) : type(2), q(q) {}
  ~BoundConds2() {}
};
// *** END OF boundary conditions (type 2) ***


// *** Boundary conditions (type 3) ***
struct BoundConds3
{
  int type;
  double T_amb, alpha;

  BoundConds3(double t_ambient_C, double alpha) :
    type(3), T_amb(t_ambient_C + T_ABS), alpha(alpha) {}
  ~BoundConds3() {}
};
// *** END OF boundary conditions (type 3) ***


#endif // TYPES_H
