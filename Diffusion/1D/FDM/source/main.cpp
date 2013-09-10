#include <initial-conditions.h>
#include <diffusion.h>
#include <utility.h>
#include <output.h>
#include <norms.h>

#include <Eigen/Dense>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace Eigen;
using namespace std;


int main() {
  const double mu      = 0.5;
  const double kappa   = 1.0;
  const double T       = 0.1;
  const int N          = 256;
  double t             = 0.0;

  #include <working-scripts/squareWave-BDF2>

  return 0;
}
