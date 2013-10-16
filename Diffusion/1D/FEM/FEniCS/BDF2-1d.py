#! /usr/bin/python

from dolfin import *
import numpy

from initialConditions import *
from utility import *

N = 256
while (N <= 8192):
  mu      = 1.0
  kappa   = 1.0
  t_end   = 0.03125

  h = 1.0 / N

  dt = mu *  h / kappa

  mesh = UnitIntervalMesh(N)
  V = FunctionSpace(mesh, 'Lagrange', 2)

  u2 = sineWave(kappa, 0)
  # u2 = squareWave()
  u1 = sineWave(kappa, 0)
  # u1 = squareWave()

  bc = DirichletBC(V, Constant(0.0), boundary)

  u_2 = interpolate(u2, V)
  u_1 = interpolate(u1, V)

  u = TrialFunction(V)
  v = TestFunction(V)

  dt_temp = dt
  a = u*v*dx + kappa * dt_temp * inner(nabla_grad(v), nabla_grad(u)) * dx
  A = assemble(a)
  u = Function(V)
  L = u_1*v*dx
  b = assemble(L)
  bc.apply(A, b)
  solve(A, u.vector(), b)
  u_1.assign(u)
  t = dt_temp

  u = TrialFunction(V)

  u = TrialFunction(V)
  a = v*u*dx + 2.0/3.0*kappa*dt*inner(nabla_grad(v), nabla_grad(u))*dx
  A = assemble(a)
  u = Function(V)

  while t < t_end:
    L = 4.0/3.0*v*u_1*dx - 1.0/3.0*v*u_2*dx
    b = assemble(L)
    bc.apply(A, b)
    solve(A, u.vector(), b)
    u_2.assign(u_1)
    u_1.assign(u)
    t += dt

  M = min(N, 300)

  u_exact = interpolate(sineWave(kappa, t), V)
  # u_exact = interpolate(fourierSquare(M, kappa, t), V)

  u_error = u_1.vector() - u_exact.vector()

  print "N: ", N, " L2: ", u_error.norm("l2")/N, " L1: ", u_error.norm("l1")/N, " Linf: ", u_error.norm("linf")
  N *= 2
