#define bc(u, N, i)   periodicBC(u, N, i)

#define min3(x, y, z)  (x < y) ? ((x < z) ? x : z) : ((y < z) ? y : z)
#define sign(x)        (x < 0) ? (-1) : 1
#define mod(x, N)      (x + N) % N

#include <Eigen/Dense>

#include <iostream>

using namespace Eigen;
using namespace std;

inline double periodicBC (Eigen::Ref<Eigen::VectorXd> u,
                          unsigned int N,
                          int i) {
  if (i < 0) return periodicBC (u, N, i + N);
  else return u[i % N];
}

inline void init_timestep (double & delta_t,
                           const double end_t) {
  int N_timestep = end_t / delta_t;
  if (abs (N_timestep * delta_t - end_t) > 10e-16) {
    delta_t = end_t / (++N_timestep);
  }
}

inline void monotonicityCheck (Eigen::VectorXd u,
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
    std::cerr << "Monotonicity not preserved!" << std::endl
              << "\tlower: " << old_lower << " -> " << new_lower << std::endl
              << "\tupper: " << old_upper << " -> " << new_upper << std::endl;
  }

  old_lower = new_lower;
  old_upper = new_upper;
}
