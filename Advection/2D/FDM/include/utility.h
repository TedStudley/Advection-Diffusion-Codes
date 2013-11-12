#define bc(u, N, i, j)   periodicBC(u, N, i, j)

#define min3(x, y, z)  (x < y) ? ((x < z) ? x : z) : ((y < z) ? y : z)
#define sign(x)        (x < 0) ? (-1) : 1
#define mod(x, N)      (x + N) % N

#define TOL 10e-6

#include <iostream>
using namespace std;

inline double periodicBC (VectorXd u,
                          const int N,
                          const int i,
                          const int j) {
  if (i < 0) return periodicBC (u, N, i + N, j);
  if (j < 0) return periodicBC (u, N, i, j + N);
  else return u[(i % N) * N + (j % N)];
}

inline void init_timestep (double & dt,
                           const double T) {
  int N_timestep = T / dt;
  if (abs (N_timestep * dt - T) > TOL) {
    dt = T / (++N_timestep);
  }
}

inline void monotonicityCheck (VectorXd u,
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
    cerr << "==> Monotonicity not preserved! <==" << endl
         << "==>\tlower: " << old_lower << " -> " << new_lower << endl
         << "==>\tupper: " << old_upper << " -> " << new_upper << endl;
  }

  old_lower = new_lower;
  old_upper = new_upper;
}
