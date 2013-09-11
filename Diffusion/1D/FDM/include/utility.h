#define bc(u, N, i)   periodicBC(u, N, i)

#define min3(x, y, z)  (x < y) ? ((x < z) ? x : z) : ((y < z) ? y : z)
#define sign(x)        (x < 0) ? (-1) : 1
#define mod(x, N)      (x + N) % N

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

inline double periodicBC (Ref<VectorXd> u,
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
