#include <norms.h>

#include <Eigen/Dense>

double maxNorm (Eigen::VectorXd u) {
  return u.array ().abs ().maxCoeff ();
}

double oneNorm (Eigen::VectorXd u) {
  return (u.array ().abs ().sum ()) / u.rows ();
}

double twoNorm (Eigen::VectorXd u) {
  return sqrt (u.array ().square ().sum () / u.rows ());
}
