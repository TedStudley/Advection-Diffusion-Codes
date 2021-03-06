#include <norms.h>

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

double maxNorm (VectorXd u) {
  return u.array ().abs ().maxCoeff ();
}

double oneNorm (VectorXd u) {
  return (u.array().abs ().sum()) / u.rows ();
}

double twoNorm (VectorXd u) {
  return sqrt (u.array ().square ().sum () / u.rows ());
}
