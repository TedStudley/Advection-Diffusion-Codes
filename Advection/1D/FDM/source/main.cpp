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

#define CHECK_MONOTONICITY true
#define dim 1

using namespace Eigen;
using namespace std;

int main() {
  int          N     = 128;
  const double mu    = 0.9;
  const double kappa = 1.0;
  const double T     = 1.0;
  Vector1d     v     = Vector1d::Constant(1.0);
  const int    k     = 2;
  double       t0    = 0.0;

  #include "workingScripts/outputScripts/sineWave-upwind"
  #include "workingScripts/outputScripts/sineWave-fromm"
  #include "workingScripts/outputScripts/sineWave-fvl"
  #include "workingScripts/outputScripts/squareWave-upwind"
  #include "workingScripts/outputScripts/squareWave-fromm"
  #include "workingScripts/outputScripts/squareWave-fvl"


  return 0;
}
