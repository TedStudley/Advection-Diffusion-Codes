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
  const double mu    = 2.0;
  const double kappa = 1.0;
  const double T     = 1.0;
  double t           = 0.0;
  const int k        = 1;
  int N              = 64;

  #include <working-scripts/OutputScripts/squareWave-bdf2> 

  return 0;
}
