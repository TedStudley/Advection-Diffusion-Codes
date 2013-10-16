def boundary(x, on_boundary):
  tol = 1E-15
  return abs(x[0]) < tol or \
         abs(x[0] - 1) < tol
