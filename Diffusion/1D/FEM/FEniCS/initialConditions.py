from dolfin import *

class squareWave(Expression):
  def eval(self, values, x):
    dx = abs(x[0] - 0.5);
    if (dx < 0.25):
      values[0] = 1
    else:
      values[0] = 0

class fourierSquare(Expression):
  def __init__(self, N, kappa, t0):
    self._N     = N
    self._t0    = t0
    self._kappa = kappa
    self._bk    = range(0, N)
    for k in range(0, N):
      self._bk[k] = (2.0 * (cos ((k + 1) * pi * 0.25) - cos ((k + 1) * pi * 0.75)) / ((k + 1) * pi) * exp (-(k + 1) * (k + 1) * self._kappa * self._t0 * pi * pi))

  def eval(self, values, x):
    values[0] = 0
    for k in range(0, self._N / 2):
      values[0] += self._bk[k] * sin ((k + 1) * pi * x[0])

class sineWave(Expression):
  def __init__(self, kappa, t0):
    self._t0    = t0
    self._kappa = kappa

  def eval(self, values, x):
    values[0] = exp (-2 * 2 * self._kappa * self._t0 * pi * pi) * sin (2 * pi * x[0])
