#! /usr/bin/python

from dolfin import *
import numpy

N       = 64
t0      = 0.0001
kappa   = 0.125
dt      = 0.001
c       = 0.5
t_end   = 1.0

mesh = UnitSquareMesh(N, N)
V = FunctionSpace(mesh, 'Lagrange', 1)

class SquareWave(Expression):
  def eval(self, values, x):
    dx = abs(x[0] - 0.5);
    dy = abs(x[1] - 0.5);
    if (dx < 0.25 and dy < 0.25):
      values[0] = 1
    else:
      values[0] = 0

class FourierSquare(Expression):
  def __init__(self, N, t0):
    self._N = N
    self._t0 = t0
    self._bk = range(0, N / 2)
    for k in range(0, N / 2):
      self._bk[k] = (2.0 * (cos ((k + 1) * pi * 0.25) - cos ((k + 1) * pi * 0.75)) / ((k + 1) * pi) * exp (-(k + 1) * (k + 1) * t0 * pi))
  def eval(self, values, x):
    values[0] = 0
    for k1 in range(0, N / 2):
      for k2 in range(0, N / 2):
        values[0] += self._bk[k1] * self._bk[k2] * sin ((k1 + 1) * pi * x[0]) * sin ((k2 + 1) * pi * x[1])

#u0 = SquareWave()
u0 = FourierSquare(N=N, t0=t0)

def boundary(x, on_boundary):
  tol = 1E-15
  return abs(x[0]) < tol or \
         abs(x[1]) < tol or \
         abs(x[0] - 1) < tol or \
         abs(x[1] - 1) < tol

bc = DirichletBC(V, Constant(0.0), boundary)

u_1 = interpolate(u0, V)

u = TrialFunction(V)
v = TestFunction(V)

a = u*v*dx+kappa*dt*0.5*inner(nabla_grad(u), nabla_grad(v))*dx

A = assemble(a)

u = Function(V)
t = dt

while t <= t_end:
  L = u_1*v*dx - kappa*dt*(1-c)*inner(nabla_grad(u_1), nabla_grad(v))*dx
  b = assemble(L)
  bc.apply(A, b)
  solve(A, u.vector(), b)
  plot(u)
  interactive()
  u_1.assign(u)
  t += dt
