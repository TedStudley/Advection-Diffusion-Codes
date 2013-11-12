#include <initial-conditions.h>
#include <advection.h>
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
  int N                = 64;
  const double mu      = 0.5;
  const double kappa   = 1.0;
  Vector2d v           = {1.0, 1.0};
  const double T       = 1.0;
  const int k1         = 2;
  const int k2         = 2;
  double t0            = 0.0;

  #include "workingScripts/outputScripts/sineWave-fvl"

  return 0;
}
