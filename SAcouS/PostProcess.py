# This file is part of PyXfem, a software distributed under the MIT license.
# For any question, please contact the authors cited below.
#
# Copyright (c) 2023
# 	Shaoqi WU <shaoqiwu@outlook.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# simple postprocessor for plotting and error computation

import numpy as np
import matplotlib.pyplot as plt

import meshio


class BasePostProcess(object):

  def __init__(self, title, *args, **kwargs):
    self.title = title

  def set_figure(self, xaxis, yaxis):
    self.fig = plt.figure()
    self.ax = self.fig.add_subplot(111)
    self.ax.set_title(self.title)
    self.ax.set_xlabel(xaxis)
    self.ax.set_ylabel(yaxis)

  def compute_error(self, sol, ana_sol, remove=None):
    # relative error

    # Compute the differences between the solutions
    differences = np.abs(ana_sol - sol)

    # Square the differences
    squared_differences = differences**2

    # Sum up the squared differences
    sum_squared_differences = np.sum(squared_differences)

    # Divide by the number of nodes
    num_nodes = ana_sol.size
    mean_squared_difference = sum_squared_differences / num_nodes

    # Compute analytical
    abs_analytical = np.sum(np.abs(ana_sol)**2) / num_nodes

    # Relative error
    l2_error = np.sqrt(mean_squared_difference / abs_analytical)

    return l2_error


class PostProcessField(BasePostProcess):

  def __init__(self, x_nodes, title, quantity='Pressure', unit='Pa'):
    super().__init__(title)
    self.x_nodes = x_nodes
    self.quantity = quantity
    self.unit = unit

  def plot_sol(self, *sols, file_name=None, save=False):
    self.set_figure('Position(m)', self.quantity + '(' + self.unit + ')')
    for sol in sols:
      self.ax.plot(self.x_nodes, sol[0], label=sol[1], linestyle=sol[2])

    self.ax.legend()
    if save:
      plt.savefig(file_name)

  def save_sol(self, *sols, file_name):
    for sol in sols:
      np.savetxt(file_name, sol[0])

  def display_layers(self, *layers_pos):
    for pos in layers_pos:
      self.ax.axvline(x=pos, ls='--', c='k')


# frequency response function postprocessor
class PostProcessFRF(BasePostProcess):

  def __init__(self, freqs, title, acoustic_indicator='SPL'):
    super().__init__(title)
    self.freqs = freqs
    self.operator = acoustic_indicator
    self.unit = ''
    if acoustic_indicator == 'SPL':
      self.unit = 'dB'

  def get_operator(self):
    if self.operator == 'SPL(dB)':
      return lambda x: 20 * np.log10(np.abs(x))
    elif self.operator == 'SPL(dB) - 2':
      print("Warning: SPL(dB) - 2 is not implemented yet!")

  def plot_sol(self, *sols, save=False, file_name=None):
    self.set_figure('Frequency(Hz)', self.operator + f'({self.unit})')
    for sol in sols:
      sol_r = self.get_operator()(sol[0])
      self.ax.plot(self.freqs, sol_r, label=sol[1], linestyle=sol[2])

    self.ax.legend()
    if save:
      plt.savefig(file_name)

  def save_sol(self, *sols, file_name):
    """
        output result as a txt file, res.txt
        {
        // header, PyacoustiX FRF result file
        // Frequency(Hz) SPL(dB)
        100 20.1
        200 30.7
        300 40.4
        400 50.1
        500 48.7
        }
        """
    with open(file_name, 'w') as file:
      # write the header
      file.write('// PyacoustiX FRF result file\n')
      file.write(f'Frequency(Hz) {self.operator}({self.unit})\n')
      for sol in sols:
        file.write(f'\n{sol[1]}')
        sol_r = self.get_operator()(sol[0])
        for i in range(len(sol_r)):
          file.write(f'{self.freqs[i]} {sol_r[i]}\n')


# write a function to plot the results on 2D/3D mesh
def plot_field(mesh, sol, title, quantity='Pressure', unit='Pa'):
  # Check if the mesh is 2D or 3D
  if mesh.dim == 2:
    # Plot the 2D mesh and the solution
    plt.figure()
    plt.title(title)
    plt.xlabel('X')
    plt.ylabel('Y')
    # import pdb
    # pdb.set_trace()
    plt.triplot(mesh.nodes.T[0], mesh.nodes.T[1], mesh.connectivity)
    plt.tricontourf(mesh.nodes.T[0], mesh.nodes.T[1], sol, cmap='jet')
    plt.colorbar(label=quantity + ' (' + unit + ')')
    plt.show()
  elif mesh.dim == 3:
    # Plot the 3D mesh and the solution
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(title)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.plot_trisurf(mesh.x_nodes,
                    mesh.y_nodes,
                    mesh.z_nodes,
                    triangles=mesh.triangles,
                    cmap='jet',
                    linewidth=0.2)
    ax.scatter(mesh.x_nodes, mesh.y_nodes, mesh.z_nodes, c=sol, cmap='jet')
    ax.colorbar(label=quantity + ' (' + unit + ')')
    plt.show()
  else:
    print("Unsupported mesh dimension.")


def save_gmsh(mesh, sol, file_name, binary):
  mesh = mesh.io_mesh
  mesh.point_data = {'Pressure': sol}
  meshio.gmsh.write(file_name, mesh, "2.2", binary)


def save_plot(mesh, sol, file_name, engine=None, binary=True):
  if engine == 'None':
    engine = 'matplotlib'
  if engine == 'gmsh':
    # write the solution to a .msh file
    save_gmsh(mesh, sol, file_name, binary)
  else:
    print("Unsupported plotting engine.")


def read_solution(file_name, read_mesh=False, engine=None):
  if engine == 'None':
    engine = 'matplotlib'
  if engine == 'gmsh':
    # read the solution from a .msh file
    mesh = meshio.read(file_name)
    if read_mesh:
      return mesh, mesh.point_data['Pressure']
    else:
      return mesh.point_data['Pressure']
  else:
    print("Unsupported plotting engine.")
