#define CHECK_MONOTONICITY true

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
  int          N     = 256;
  const double mu    = 0.75;
  const double kappa = 1.0;
  const double T     = 1.0;
  Vector1d     v     = Vector1d::Constant(1.0);
  const int    k     = 2;
  double       t0    = 0.0;

  #include "workingScripts/latexScripts/sineWave-beamWarming"

  return 0;
}
