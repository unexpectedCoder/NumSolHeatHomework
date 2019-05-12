#ifndef TYPES_H
#define TYPES_H

#include <fstream>
#include <vector>
#include <string>
#include <math.h>

#define EPS 1e-12       // Accuracy
#define T_ABS 273.15    // Absolute difference between C and K


const double g = 9.80665;
const double C = 5.67;


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

  size_t N;               // Number of spacing segments
  double r1, r2;          // Inner and outer radius of cylinder
  double step;            // Space step
  double *T_table;        // Table temperature
  double *lambda;         // lambda(T) - wall's heat transfer coeff
  size_t dataSize;        // Size of lambda data
  double *T;              // T(r): [0][:] - coordinates, [1][:] - temperature
  double *r;              // Coordinates
  double rho;             // Material density (T)
  double *c;              // Specific heat (T)

  bool is_lambda;         // Was the lambda(T) initialized
  bool is_T;              // Was the T(tau) initialized
  bool is_rho;            // Was the rho(T) initialized
  bool is_c;              // Was the c(T) initialized
  bool is_grid;           // Was the r initialized
  bool is_T_table;        // Was the table temperature initialized

  double epsilon;         // Blackness
  std::string material;   // Material's name

  Wall(double r1, double r2, size_t n, const std::string &material = "NONE");
  Wall(const Wall &w);
  ~Wall();

  Wall& operator=(const Wall &w);

  void setGrid();
  void setLambdaT(const std::string& file_path);
  void setLambdaT(const double *T, const double *lam, size_t n);
  void setBlackness(double epsilon);
  void setDens(double rho);
  void setSpecificHeat(double *c, size_t n);

  inline friend std::ostream& operator<<(std::ostream &os, const Wall &w);
};


std::ostream& operator<<(std::ostream &os, const Wall &w)
{
  os << "\n\t..... WALL .....\n"
     << "\tmaterial: " << w.material << '\n'
     << "\tN: " << w.N << '\n'
     << "\tr1, m: " << w.r1 << '\n'
     << "\tr2, m: " << w.r2 << '\n'
     << "\tstep (dr), m: " << w.step << '\n'
     << "\tepsilon: " << w.epsilon << '\n'
     << "\trho, kg/m^3: " << w.rho << '\n';

  if (w.is_c)
  {
    os << "\tt, C\tc\n";
    for (size_t i = 0; i < w.dataSize; ++i)
      os << '#' << i + 1 << ".\t"
         << w.T_table[i] - T_ABS << '\t' << w.c[i] << '\n';
  }

  if (w.is_T_table)
  {
    os << "\tlambda(T):\n"
       << "\t\tT, C\tlambda, W/(m*K)\n";
    for (size_t i = 0; i < w.dataSize; ++i)
      os << "\t#" << i + 1 << ".\t"
         << w.T_table[i] - T_ABS << '\t' << w.lambda[i] << '\n';
  }

  if (w.is_T)
  {
    os << "\tT(r):\n\tr, m\tT, C\n";
    for (size_t i = 0; i < w.N; ++i)
      os << '#' << i + 1 << ".\t"
         << w.r[i] << '\t' << w.T[i] - T_ABS << '\n';
  }

  os << '\n';
  return os;
}


typedef std::vector<Wall> Walls;
typedef std::vector<Wall>::iterator WallItr;
typedef std::vector<Wall>::const_iterator WallCItr;
// *** END OF Wall ***


// *** Environment ***
struct Environment
{
  double Ta;        // Ambient temperature
  // Table data for future interpolation
  size_t dataSize;
  double *T;
  double *lambda;
  double *a;
  double *c;
  double *rho;
  double *nu;
  double *mu;
  double *Pr;

  Environment() : Ta(0.0), dataSize(0) {}
  ~Environment();

  inline friend std::ostream& operator<<(std::ostream &os,
                                         const Environment &e);
};


std::ostream& operator<<(std::ostream &os, const Environment &e)
{
  os << "\t..... ENVIRONMENT .....\n";
  os << "\tt, C\t\tlambda\t\trho\t\tcp\t\ta\t\tnu\t\tmu\t\tPr\n";
  for (size_t i = 0; i < e.dataSize; ++i)
    os << '\t' << e.T[i] - T_ABS
       << "\t\t" << e.lambda[i]
       << "\t\t" << e.rho[i]
       << "\t\t" << e.c[i]
       << "\t\t" << e.a[i]
       << "\t\t" << e.nu[i]
       << "\t\t" << e.mu[i]
       << "\t\t" << e.Pr[i] << '\n';

  os << '\n';
  return os;
}
// *** END OF Environment ***


// *** Boundary conditions (type 2) ***
struct BoundCond
{
  int type;             // Type of boundary conditions:
  double q;             // for type 2 (heat flow to wall)
  double T_amb, alpha;  // for type 3 (ambient T and heat emission coeff)

  BoundCond() : type(0), q(0.0), T_amb(0.0), alpha(0.0) {}

  void setType2(double q)
  {
    type = 2;
    this->q = q;
  }
  void setType3(double t_amb_C, double alpha = 0.0)
  {
    type = 3;
    T_amb = t_amb_C + T_ABS;
    this->alpha = alpha;
  }
};
// *** END OF boundary conditions ***


// *** Starting conditions ***
struct StartConds
{
  Error err;

  double T0;    // Walls starting temperature
  double D, H;  // Cylinder geometry
  double time;  // Starting time

  StartConds(double t_C, double t = 0.0) :
    T0(t_C + T_ABS), time(t) {}
  ~StartConds() {}

  void setGeometry(const Walls &ws, double h);
};
// *** END OF Starting conditions ***


#endif // TYPES_H
