#define bc(u, N, i)   periodicBC(u, N, i)

#define min3(x, y, z)  (x < y) ? ((x < z) ? x : z) : ((y < z) ? y : z)
#define sign(x)        (x < 0) ? (-1) : 1
#define mod(x, N)      (x + N) % N

#include <iostream>
using namespace std;

inline double periodicBC (Ref<VectorXd> u,
                          unsigned int N,
                          int i) {
  if (i < 0) return periodicBC (u, N, i + N);
  else return u[i % N];
}
