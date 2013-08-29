#define bc(u, N, i)   periodicBC(u, N, i)

#define min3(x, y, z)  (x < y) ? ((x < z) ? x : z) : ((y < z) ? y : z)
#define sign(x)        (x < 0) ? (-1) : ((x == 0) ? 0 : 1)

inline double periodicBC (Ref<VectorXd> u,
                          unsigned int N,
                          int i) {
  while (i < 0) i += N;
  while (i >= int(N)) i -= N;
  return u[i];
}

inline double diricheletBC (Ref<VectorXd> u,
                            unsigned int N,
                            int i) {
  if (i <= 0) return u[0];
  if (i >= int (N - 1)) return u[N-1];
  return u[i];
}
