#include <norms.h>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

double maxNorm (VectorXd u) {
  return u.array ().abs ().maxCoeff ();
}

double oneNorm (VectorXd u) {
  return (u.array ().abs ().sum ()) / u.rows ();
}

double twoNorm (VectorXd u) {
  return sqrt (u.array ().square ().sum () / u.rows ());
}
