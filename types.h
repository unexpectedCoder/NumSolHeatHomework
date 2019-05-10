#ifndef TYPES_H
#define TYPES_H

#include <fstream>
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


// *** Wall (settlement layer) ***
struct Wall
{
  Error err;

  size_t N;                  // Number of spacing segments
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

  Wall(double r1, double r2, size_t n, const std::string &material = "NONE");
  Wall(const Wall &w);
  ~Wall();

  Wall& operator=(const Wall &w);

  void setLambda(const std::string& file_path);
  void setLambda(const double *T, const double *lam, size_t n);
  void setStartTemperature(double temp_C);
  void setBlackness(double epsilon);
  void setDens(double rho);
  void setSpecificHeat(double c);

  inline friend std::ostream& operator<<(std::ostream &os, const Wall &w);
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

  os << "\tT(r):\n\tr, m\tT, C\n";
  if (w.T)
    for (size_t i = 0; i < w.N; ++i)
      if (w.T[i])
        os << '\t' << w.T[0][i] << '\t' << w.T[1][i] - T_ABS << '\n';

  os << '\n';
  return os;
}


typedef std::vector<Wall> Walls;
typedef std::vector<Wall>::iterator WallItr;
typedef std::vector<Wall>::const_iterator WallCItr;
// *** END OF Wall ***


// *** Boundary conditions (type 2) ***
struct BoundConds
{
  int type;             // Type of boundary conditions:
  double q;             // for type 2 (heat flow to wall)
  double T_amb, alpha;  // for type 3 (ambient T and heat emission coeff)

  BoundConds() : type(0), q(0.0), T_amb(0.0), alpha(0.0) {}
};
// *** END OF boundary conditions ***


// *** Starting conditions ***
struct StartConds
{
  double T, T_amb;
  double D, H;
  double time;

  StartConds(double t_C, double t_amb_C, double d, double h, double t = 0.0) :
    T(t_C + T_ABS), T_amb(t_amb_C + T_ABS), D(d), H(h), time(t) {}
  ~StartConds() {}
};
// *** END OF Starting conditions ***


#endif // TYPES_H
