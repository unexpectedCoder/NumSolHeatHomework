#include "types.h"

#include <iostream>

using namespace std;


// *** Wall ***
Wall::Wall(double r1, double r2, size_t n, const std::string &material) :
  dataSize(0),
  is_lambda(false), is_T(false), is_rho(false),
  is_c(false), is_r(false), is_T_table(false),
  epsilon(1.0), material(material)
{
  if (r1 < 0.0 || fabs(r2 - r1) < EPS)
    throw err.sendEx("r1 < 0 or r2 < r1");
  if (n < 2)
    throw err.sendEx("number of segments < 1");

  this->r1 = r1;
  this->r2 = r2;
  N = n + 1;
  step = (r2 - r1) / (N - 1);

  r = new double[N];
  for (size_t i = 0; i < N; ++i)
    r[i] = step * i;
  is_r = true;
  // This array is used in vector, so memory is fully allocated here
  T = new double[N];
}


Wall::Wall(const Wall &w)
{
  /*
   * Copy constructor is needed for the correct
   * deliting ** pointers of this structure
   * in case of using std::vector
   * or using this structure in other parts of code.
  */

  N = w.N;
  r1 = w.r1;
  r2 = w.r2;
  step = w.step;
  dataSize = w.dataSize;
  epsilon = w.epsilon;
  material = w.material;
  rho = w.rho;

  is_lambda = w.is_lambda;
  is_rho = w.is_rho;
  is_c = w.is_c;
  is_T = w.is_T;
  is_r = w.is_r;

  T_table = new double[dataSize];
  lambda = new double[dataSize];
  c = new double[dataSize];
  for (size_t i = 0; i < dataSize; ++i)
  {
    T_table[i] = w.T_table[i];
    lambda[i] = w.lambda[i];
    c[i] = w.c[i];
  }

  r = new double[N];
  T = new double[N];
  if (is_T)
    for (size_t i = 0; i < N; ++i)
    {
      r[i] = w.r[i];
      T[i] = w.T[i];
    }
}


Wall::~Wall()
{
  delete [] lambda;
  delete [] c;
  delete [] T_table;

  if (is_r)
    delete [] r;
  if (is_T)
    delete [] T;
}


Wall& Wall::operator=(const Wall &w)
{
  /*
   * Overwriting this operator is needed
   * for the same reasons as the copy constructor.
  */

  if (this == &w)
    return *this;

  N = w.N;
  r1 = w.r1;
  r2 = w.r2;
  step = w.step;
  dataSize = w.dataSize;
  epsilon = w.epsilon;
  material = w.material;
  rho = w.rho;

  is_lambda = w.is_lambda;
  is_rho = w.is_rho;
  is_c = w.is_c;
  is_T = w.is_T;
  is_r = w.is_r;

  T_table = new double[dataSize];
  lambda = new double[dataSize];
  c = new double[dataSize];
  for (size_t i = 0; i < dataSize; ++i)
  {
    T_table[i] = w.T_table[i];
    lambda[i] = w.lambda[i];
    c[i] = w.c[i];
  }

  r = new double[N];
  T = new double[N];
  if (is_T)
    for (size_t i = 0; i < N; ++i)
    {
      r[i] = w.r[i];
      T[i] = w.T[i];
    }

  return *this;
}


void Wall::setLambdaT(const string& file_path)
{
  /*
   * Read function lambda(T) from source file
   * and set table **<var>[2][size] with T and lambda.
  */

  if (is_lambda)
    throw err.sendEx("lamba is already set!");
  is_lambda = true;

  ifstream file(file_path.c_str(), ios_base::in);
  if (!file.is_open())
    throw err.sendEx("file is not opened");

  double buf;
  while (!file.eof())
  {
    file >> buf >> buf;
    dataSize++;
  }
  dataSize--;
  file.close();

  file.open(file_path.c_str(), ios_base::in);
  lambda = new double[dataSize];
  T_table = new double[dataSize];
  for (size_t i = 0; i < dataSize; ++i)
  {
    file >> T_table[i] >> lambda[i];
    T_table[i] += T_ABS;
  }
  file.close();

  is_T_table = true;
}


void Wall::setLambdaT(const double *T, const double *lam, size_t n)
{
  /*
   * Set table **<var>[2][size] with T and lambda
   * using source array T and lambda.
  */

  if (dataSize != 0 && n != dataSize)
    throw err.sendEx("set size != data size");
  if (n < 2)
    throw err.sendEx("arrays' size < 2");
  if (is_lambda)
    throw err.sendEx("lamba is already set");
  is_lambda = true;

  if (dataSize == 0)
    dataSize = n;

  lambda = new double[dataSize];
  T_table = new double[dataSize];
  for (size_t i = 0; i < dataSize; ++i)
  {
    T_table[i] = T[i] + T_ABS;
    lambda[i] = lam[i];
  }

  is_T_table = true;
}


void Wall::setBlackness(double epsilon)
{
  if (epsilon < 0.0 || epsilon > 1.0)
    throw err.sendEx("blackness must be in [0; 1]");
  this->epsilon = epsilon;
}


void Wall::setDens(double rho_)
{
  if (rho_ < 0.0)
    throw err.sendEx("density must be greater than 0");
  rho = rho_;
  is_rho = true;
}


void Wall::setSpecificHeat(double *c_, size_t n)
{
  if (n != dataSize && dataSize != 0)
    throw err.sendEx("size of c != size of data");
  for (size_t i = 0; i < n; ++i)
    if (c_[i] < 0.0)
      throw err.sendEx("specific heat must be greater than 0");

  c = new double[dataSize];
  for (size_t i = 0; i < dataSize; ++i)
    c[i] = c_[i];
  is_c = true;
}
// *** END OF Wall ***


// *** Environment ***
Environment::~Environment()
{
  if (dataSize != 0)
  {
    delete [] T;
    delete [] lambda;
    delete [] a;
    delete [] c;
    delete [] rho;
    delete [] nu;
    delete [] mu;
    delete [] Pr;
  }
}
// *** END OF Environment ***


// *** Start conds ***
void StartConds::setGeometry(const Walls &ws, double h)
{
  size_t n = ws.size();
  H = h;

  if (n < 2)
  {
    D = 2.0 * ws[0].r2;
    return;
  }

  for (size_t i = 0; i < n - 1; ++i)
    if (fabs(ws[i + 1].r1 - ws[i].r2) > EPS)
      throw err.sendEx("cylinder diameter is not equal to walls sizes");
  D = 2.0 * ws[n - 1].r2;
}
// *** END OF Start conds ***
