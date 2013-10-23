#define bc(u, N, i)   periodicBC(u, N, i)

#define min3(x, y, z)  (x < y) ? ((x < z) ? x : z) : ((y < z) ? y : z)
#define sign(x)        (x < 0) ? (-1) : 1
#define mod(x, N)      (x + N) % N

#define TOL 10e-6

#include <iostream>

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

typedef Matrix<double, 1, 1> Vector1d;

inline double periodicBC (Ref<Eigen::VectorXd> u,
                          unsigned int         N,
                          int                  i) {
  if (i < 0) return periodicBC (u, N, i + N);
  else return u[i % N];
}

inline void init_timestep (double &     dt,
                           const double T) {
  int N_timestep = T / dt;
  if (abs (N_timestep * dt - T) > TOL) {
    dt = T / (++N_timestep);
  }
}

inline void monotonicityCheck (VectorXd     u,
                               const double lower = 0.0,
                               const double upper = 1.0) {
  static double old_lower;
  static double old_upper;

  static bool next_step;
  if (!next_step) {
    old_lower = lower;
    old_upper = upper;
    next_step = true;
  }

  double new_lower = u.minCoeff ();
  double new_upper = u.maxCoeff ();

  if (new_lower < old_lower || new_upper > old_upper) {
    std::cerr << "==> Monotonicity not preserved! <==" << endl
              << "==>\tlower: " << old_lower << " -> " << new_lower << endl
              << "==>\tupper: " << old_upper << " -> " << new_upper << endl;
  }

  old_lower = new_lower;
  old_upper = new_upper;
}
