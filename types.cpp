#include "types.h"


using namespace std;


// *** Wall ***
Wall::Wall(double r1, double r2, size_t n, const std::string &material) :
  dataSize(0), is_lambda(false), is_crho(false), is_T(false),
  epsilon(1.0), rho(0.0), c(0.0), material(material)
{
  if (r1 < 0.0 || fabs(r2 - r1) < EPS)
    throw err.sendEx("r1 < 0 or r2 < r1");
  if (n < 2)
    throw err.sendEx("number of segments < 1");

  this->r1 = r1;
  this->r2 = r2;
  N = n + 1;
  step = (r2 - r1) / (N - 1);

  lambda_T = new double*[2];
  crho = new double*[2];

  // Especially selected T, since filled inside the vector
  T = new double*[2];
  for (size_t i = 0; i < 2; ++i)
    T[i] = new double[N];
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
  rho = w.rho;
  c = w.c;
  material = w.material;

  is_lambda = w.is_lambda;
  is_crho = w.is_crho;
  is_T = w.is_T;

  lambda_T = new double*[2];
  crho = new double*[2];
  for (size_t i = 0; i < 2; ++i)
  {
    lambda_T[i] = new double[dataSize];
    crho[i] = new double[dataSize];
  }
  for (size_t i = 0; i < dataSize; ++i)
  {
    lambda_T[0][i] = w.lambda_T[0][i];
    lambda_T[1][i] = w.lambda_T[1][i];
    crho[0][i] = w.crho[0][i];
    crho[1][i] = w.crho[1][i];
  }

  // Especially selected T, since filled inside the vector
  T = new double*[2];
  for (size_t i = 0; i < 2; ++i)
    T[i] = new double[N];
  if (is_T)
    for (size_t i = 0; i < N; ++i)
    {
      T[0][i] = w.T[0][i];
      T[1][i] = w.T[1][i];
    }
}


Wall::~Wall()
{
  if (is_lambda)
    for (size_t i = 0; i < 2; ++i)
      delete [] lambda_T[i];
  if (is_crho)
    for (size_t i = 0; i < 2; ++i)
      delete [] crho[i];
  delete [] lambda_T;
  delete [] crho;

  for (size_t i = 0; i < 2; ++i)
    delete [] T[i];
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
  rho = w.rho;
  c = w.c;
  material = w.material;

  is_lambda = w.is_lambda;
  is_crho = w.is_crho;
  is_T = w.is_T;

  lambda_T = new double*[2];
  crho = new double*[2];
  for (size_t i = 0; i < 2; ++i)
  {
    lambda_T[i] = new double[dataSize];
    crho[i] = new double[dataSize];
  }

  for (size_t i = 0; i < dataSize; ++i)
  {
    lambda_T[0][i] = w.lambda_T[0][i];
    lambda_T[1][i] = w.lambda_T[1][i];
    crho[0][i] = w.crho[0][i];
    crho[1][i] = w.crho[1][i];
  }

  // Especially selected T, since filled inside the vector
  T = new double*[2];
  for (size_t i = 0; i < 2; ++i)
    T[i] = new double[N];
  if (is_T)
    for (size_t i = 0; i < N; ++i)
    {
      T[0][i] = w.T[0][i];
      T[1][i] = w.T[1][i];
    }

  return *this;
}


void Wall::setLambda(const string& file_path)
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
    ++dataSize;
  }
  dataSize--;
  file.close();

  file.open(file_path.c_str(), ios_base::in);
  for (size_t i = 0; i < 2; ++i)
    lambda_T[i] = new double[dataSize];
  for (size_t i = 0; i < dataSize; ++i)
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

  if (dataSize != 0 && n != dataSize)
    throw err.sendEx("set size != data size");
  if (n < 2)
    throw err.sendEx("arrays' size < 2");
  if (is_lambda)
    throw err.sendEx("lamba is already set");
  is_lambda = true;

  if (dataSize == 0)
    dataSize = n;
  for (size_t i = 0; i < 2; ++i)
    lambda_T[i] = new double[dataSize];
  for (size_t i = 0; i < dataSize; ++i)
  {
    lambda_T[0][i] = T[i];
    lambda_T[1][i] = lam[i];
  }
}


void Wall::set_crho(const double *T, const double *c_rho, size_t n)
{
  if (dataSize != 0 && n != dataSize)
    throw err.sendEx("set size != data size");
  if (n < 2)
    throw err.sendEx("arrays' size < 2");
  if (is_crho)
    throw err.sendEx("c*rho is already set");
  is_crho = true;

  if (dataSize == 0)
    dataSize = n;
  for (size_t i = 0; i < 2; ++i)
    crho[i] = new double[dataSize];
  for (size_t i = 0; i < dataSize; ++i)
  {
    crho[0][i] = T[i];
    crho[1][i] = c_rho[i];
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
