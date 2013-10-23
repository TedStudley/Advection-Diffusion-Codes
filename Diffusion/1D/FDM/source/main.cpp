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

#define CHECK_MONOTONICITY true
#define DIM 1

using namespace Eigen;
using namespace std;


int main() {
  int          N       = 256;
  const double mu      = 0.9;
  const double kappa   = 1.0;
  const double T       = 0.25;
  Vector1d     v       = Vector1d::Constant(1.0);
  const int k          = 2;
  double t0            = 0.0;

  #include <workingScripts/latexScripts/sineWave-BDF2>

  return 0;
}
