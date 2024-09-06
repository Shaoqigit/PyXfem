import numpy as np
from abc import ABCMeta, abstractmethod
from cmath import sin, cos


class ImpedenceKundltTube():

  def __init__(self, *args) -> None:
    self.mat = args[0]
    self.omega = args[1]
    self.nature_bc = args[2]
    assert (self.nature_bc['type'] == 'fluid_velocity')
    self.impedence_bc = args[3]
    assert (self.impedence_bc['type'] == 'impedence')
    self.analytical_field()

  def analytical_field(self):
    Z_0 = self.mat.Z_f
    Up = self.nature_bc['value']
    A = self.impedence_bc['value']
    L = 1
    k = self.omega / self.mat.c_f

    alpha = Z_0 * Up * (1 + A) / ((1 + A) * np.exp(1j * k * L) -
                                  (1 - A) * np.exp(-1j * k * L))
    beta = alpha * (1 - A) / (1 + A)
    self.p = lambda x: alpha * np.exp(-1j * k * x) + beta * np.exp(1j * k * x)

  def sol_on_mesh(self, mesh, sol_type='pressure'):
    # import pdb; pdb.set_trace()
    ana_sol = np.zeros(
        mesh.nb_nodes,
        dtype=np.complex128)    #initialize the analytical solution vector
    for i, x in enumerate(mesh.nodes):
      ana_sol[i] = self.p(x)
    return ana_sol


class DoubleleLayerKundltTube():

  def __init__(self, *args) -> None:
    self.l1 = args[0]
    self.l2 = args[1]
    self.mat1 = args[2]
    self.mat2 = args[3]
    self.omega = args[4]
    self.bc = args[5]
    assert (self.bc['type'] == 'fluid_velocity')
    self.analytical_field()

  def analytical_field(self):
    l_a = self.l1
    l_eq = self.l2
    v_0 = self.bc['value']
    k_a = self.omega / self.mat1.c_f
    rho_a = self.mat1.rho_f

    self.mat2.set_frequency(self.omega)
    k_eq = self.omega / self.mat2.c_f
    rho_eq = self.mat2.rho_f

    Thickness, sigma_film = 0, 0

    self.P_analy = []
    self.v_analy = []
    C = complex(-k_a**2 * k_eq * sin(k_a * l_a) * sin(k_eq * l_eq) *
                sigma_film * Thickness + 1j * self.omega * k_a**2 * rho_eq *
                sin(k_a * l_a) * cos(k_eq * l_eq) + 1j * self.omega * k_a *
                k_eq * rho_a * sin(k_eq * l_eq) * cos(k_a * l_a))

    self.P_analy.append(lambda *x: (
        (1j * v_0 * sigma_film * Thickness * self.omega * k_a * k_eq * rho_a *
         sin(k_eq * l_eq) + v_0 * self.omega**2 * rho_a * rho_eq * k_a * cos(
             k_eq * l_eq)) * cos(k_a * x[0]) / C + v_0 * self.omega**2 * rho_a
        **2 * k_eq * sin(k_eq * l_eq) * sin(k_a * x[0]) / C))
    self.v_analy.append(lambda *x: (
        (1j * v_0 * sigma_film * Thickness * self.omega * k_a * k_eq * rho_a *
         sin(k_eq * l_eq) + v_0 * self.omega**2 * rho_a * rho_eq * k_a * cos(
             k_eq * l_eq)) * k_a * sin(k_a * x[0]) /
        (1j * self.omega * rho_a * C) + v_0 * self.omega**2 * rho_a**2 * k_eq *
        sin(k_eq * l_eq) * k_a * cos(k_a * x[0]) /
        (-1j * self.omega * rho_a * C)))

    self.P_analy.append(
        lambda *x: (v_0 * self.omega**2 * rho_a * rho_eq * k_a * cos(
            k_eq * l_eq) * cos(k_eq * x[0]) / C + v_0 * self.omega**2 * rho_a *
                    rho_eq * k_a * sin(k_eq * l_eq) * sin(k_eq * x[0]) / C))
    self.v_analy.append(lambda *x: (
        v_0 * self.omega**2 * rho_a * rho_eq * k_a * cos(k_eq * l_eq) * k_eq *
        sin(k_eq * x[0]) / (1j * self.omega * rho_a * C) + v_0 * self.omega**2
        * rho_a * rho_eq * k_a * sin(k_eq * l_eq) * k_eq * cos(k_eq * x[0]) /
        (-1j * self.omega * rho_a * C)))

  def sol_on_mesh(self, mesh, sol_type):
    ana_sol = np.zeros(mesh.nb_nodes, dtype=np.complex128)
    if mesh.dim == 1:
      for i, x in enumerate(mesh.nodes):
        if x <= 0:
          if sol_type == 'pressure':
            sol = self.P_analy[0](x)
          elif sol_type == 'fluid_velocity':
            sol = self.v_analy[0](x)
        elif x >= 0:
          if sol_type == 'pressure':
            sol = self.P_analy[1](x)
          elif sol_type == 'fluid_velocity':
            sol = self.v_analy[1](x)
        ana_sol[i] = sol

    elif mesh.dim >= 2:
      for i, x in enumerate(mesh.nodes):
        if x[0] <= 0:
          if sol_type == 'pressure':
            sol = self.P_analy[0](*x)
          elif sol_type == 'fluid_velocity':
            sol = self.v_analy[0](*x)
        elif x[0] >= 0:
          if sol_type == 'pressure':
            sol = self.P_analy[1](*x)
          elif sol_type == 'fluid_velocity':
            sol = self.v_analy[1](*x)
        ana_sol[i] = sol
    return ana_sol


class ObliquePlaneWave():

  def __init__(self, *args) -> None:
    self.mat1 = args[0]
    self.mat2 = args[1]
    self.omega = args[2]
    self.angle = args[3]
    self.A = args[4]
    self.analytical_field()

  def analytical_field(self):
    Z_0 = self.mat1.Z_f
    self.mat2.set_frequency(self.omega)
    k1 = self.omega / self.mat1.c_f
    k2 = self.omega / self.mat2.c_f
    # convert to radian
    theta = self.angle * np.pi / 180
    k1_x = k1 * np.cos(theta)
    k1_y = k1 * np.sin(theta)
    k2_x = np.sqrt(pow(k2, 2) - pow(k1, 2) * pow(np.sin(theta), 2))
    k2_y = k1_y
    rho_1 = self.mat1.rho_f
    rho_2 = self.mat2.rho_f
    T = 2. / (1. + rho_1 * k2_x / (rho_2 * k1_x))
    R = 1. - rho_1 * k2_x / (rho_2 * k1_x) * T

    self.p = []
    self.p.append(lambda x, y: np.exp(-1j * k1_x * x - 1j * k1_y * y) + R * np.
                  exp((1j * k1_x * x - 1j * k1_y * y)))
    self.p.append(lambda x, y: T * np.exp(-1j * k2_x * x - 1j * k2_y * y))

    self.v = []    #v_x(<0), v_x(>0), v_y(<0), v_y(>0)
    self.v.append(
        lambda x, y: -1j * k1_x * np.exp(-1j * k1_x * x - 1j * k1_y * y) + 1j *
        k1_x * R * np.exp(1j * k1_x * x - 1j * k1_y * y))
    self.v.append(
        lambda x, y: -1j * k2_x * T * np.exp(-1j * k2_x * x - 1j * k2_y * y))
    self.v.append(
        lambda x, y: -1j * k1_y * np.exp(-1j * k1_x * x - 1j * k1_y * y) + 1j *
        k1_y * R * np.exp(1j * k1_x * x - 1j * k1_y * y))
    self.v.append(
        lambda x, y: -1j * k2_y * T * np.exp(-1j * k2_x * x - 1j * k2_y * y))

    # x<0
    self.velocity_1 = lambda x, y: np.array([
        1 / rho_1 *
        (-1j * k1_x * np.exp(-1j * k1_x * x - 1j * k1_y * y) + 1j * k1_x * R *
         np.exp(1j * k1_x * x - 1j * k1_y * y)), 1 / rho_1 *
        (-1j * k1_y * np.exp(-1j * k1_x * x - 1j * k1_y * y) - 1j * k1_y * R *
         np.exp(1j * k1_x * x - 1j * k1_y * y))
    ])
    self.velocity_2 = lambda x, y: np.array([
        1 / rho_2 * -1j * k2_x * T * np.exp(-1j * k2_x * x - 1j * k2_y * y), 1
        / rho_2 * -1j * k2_y * T * np.exp(-1j * k2_x * x - 1j * k2_y * y)
    ])
    # self.velocity_2 = lambda x, y: np.array([
    #     1 / rho_2 * -1j * k2_x * T * np.exp(-1j * k2_x * x - 1j * k2_y * y), 0.
    # ])

  def sol_on_mesh(self, mesh, sol_type='pressure'):
    ana_sol = np.zeros((mesh.nb_nodes), dtype=np.complex128)
    for i, x in enumerate(mesh.nodes):
      if x[0] <= 0:
        if sol_type == 'pressure':
          sol = self.p[0](x[0], x[1])
        elif sol_type == 'fluid_velocity':
          sol = (self.v[0](x[0], x[1]), self.v[2](x[0], x[1]))
      elif x[0] >= 0:
        if sol_type == 'pressure':
          sol = self.p[1](x[0], x[1])
        elif sol_type == 'fluid_velocity':
          sol = (self.v[1](x[0], x[1]), self.v[3](x[0], x[1]))
      ana_sol[i] = sol
    return ana_sol


class MonopoleSourceSolution():

  def __init__(self, mesh, *args) -> None:
    self.mesh = mesh
    self.mat = args[0]
    self.omega = args[1]
    self.bc = args[2]
    assert (self.bc['type'] == 'fluid_velocity')
    self.analytical_field()

  def analytical_field(self):
    Z_0 = self.mat.Z_f
    Up = self.bc['value']
    k = self.omega / self.mat.c_f
    rho = self.mat.rho_f
    self.p = lambda x: Up * np.exp(-1j * k * x) / (1j * k * rho)

  def sol_on_nodes(self, ana_sol, sol_type='pressure'):
    for i, x in enumerate(self.mesh.nodes):
      ana_sol[i] = self.p(x)
    return ana_sol


if __name__ == '__main__':
  import numpy as np
  import matplotlib.pyplot as plt
  from scipy.special import hankel1

  def gaussian_source(x, y, Sx, Sy, alpha):
    return (1 / (np.abs(alpha) * np.sqrt(np.pi))) * np.exp(-(((x - Sx)**2 +
                                                              (y - Sy)**2) /
                                                             (alpha**2)))

  def green_function(x, y, x_prime, y_prime, k):
    r = np.sqrt((x - x_prime)**2 + (y - y_prime)**2)
    return (1j / 4) * hankel1(0, k * r)

  def helmholtz_solution(x, y, Sx, Sy, alpha, k, grid_size=100):
    x_prime = np.linspace(-0.5, 1.5, grid_size)
    y_prime = np.linspace(-0.5, 1.5, grid_size)
    X_prime, Y_prime = np.meshgrid(x_prime, y_prime)

    f = gaussian_source(X_prime, Y_prime, Sx, Sy, alpha)
    u = np.zeros_like(x, dtype=np.complex128)

    dx = x_prime[1] - x_prime[0]
    dy = y_prime[1] - y_prime[0]

    for i in range(len(x)):
      for j in range(len(y)):
        integrand = green_function(x[i], y[j], X_prime, Y_prime, k) * f
        u[i, j] = np.sum(integrand) * dx * dy

    return np.real(u)

  # Parameters
  Sx, Sy = 0.5, 0.5    # Source position
  alpha = 0.001    # Width of the Gaussian source
  k = 1    # Wave number
  grid_size = 100    # Resolution of the grid

  # Field points
  x = np.linspace(0, 1, grid_size)
  y = np.linspace(0, 1, grid_size)
  X, Y = np.meshgrid(x, y)

  # Calculate solution
  u = helmholtz_solution(X, Y, Sx, Sy, alpha, k, grid_size)

  # Plot the solution
  plt.figure(figsize=(8, 6))
  plt.contourf(X, Y, u, levels=50, cmap='jet')
  plt.colorbar(label='Acoustic Pressure')
  plt.xlabel('x')
  plt.ylabel('y')
  plt.title('Analytical Solution of Helmholtz Equation')
  plt.show()
