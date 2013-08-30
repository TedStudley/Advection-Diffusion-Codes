#! /usr/bin/python

from dolfin import *
import numpy

N       = 256
t0      = 0.0001
kappa   = 0.125
dt      = 0.001
c       = 0.5
t_end   = 1.0

mesh = UnitIntervalMesh(N)
V = FunctionSpace(mesh, 'Lagrange', 1)

class SquareWave(Expression):
  def eval(self, values, x):
    dx = abs(x[0] - 0.5);
    if (dx < 0.25):
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
    for k in range(0, N / 2):
      values[0] += self._bk[k] * sin ((k + 1) * pi * x[0])

u0 = SquareWave()
#u0 = FourierSquare(N=N, t0=t0)

def boundary(x, on_boundary):
  tol = 1E-15
  return abs(x[0]) < tol or \
         abs(x[0] - 1) < tol

bc = DirichletBC(V, Constant(0.0), boundary)

print "==> Interpolating initial data"

u_1 = interpolate(u0, V)

print "==> Setting up trial and test functions."

u = TrialFunction(V)
v = TestFunction(V)

print "==> discretizing left-hand side"

a = u*v*dx+kappa*dt*0.5*inner(nabla_grad(u), nabla_grad(v))*dx

print "==> Assembling left-hand side"

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
