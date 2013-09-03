#define bc(u, N, i, j)   periodicBC(u, N, i, j)

#define min3(x, y, z)  (x < y) ? ((x < z) ? x : z) : ((y < z) ? y : z)
#define sign(x)        (x < 0) ? (-1) : 1
#define mod(x, N)      (x + N) % N

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
