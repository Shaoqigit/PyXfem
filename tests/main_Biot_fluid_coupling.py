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

# Main Test case
import sys

sys.path.append('/home/shaoqi/Devlop/PyXfem/PyAcoustiX/')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import spy

from SAcouS.acxfem import Lobbato1DElement
from SAcouS.Mesh import Mesh1D
from SAcouS.acxfem import DofHandler1D, GeneralDofHandler1D, FESpace
from SAcouS.acxfem import HelmholtzAssembler, BiotAssembler, CouplingAssember
from SAcouS.acxfem import ApplyBoundaryConditions
from SAcouS.Materials import Air, Fluid, EquivalentFluid, PoroElasticMaterial
from SAcouS.acxfem import check_material_compability, display_matrix_in_array, plot_matrix_partten
from SAcouS.acxfem import LinearSolver
from SAcouS.acxfem import PostProcessField
from analytical.Fluid_Biot_sol import Fluid_Biot_Pressure, Fluid_Biot_Displacement, u_a, u_t, sigma_xx, P_a


def test_case():
  num_elem = 400    # number of elements
  num_nodes = num_elem + 1    # number of nodes

  nodes = np.linspace(-1, 1, num_nodes)
  # print()
  elem_connec1 = np.arange(0, num_elem)
  elem_connec2 = np.arange(1, num_nodes)
  connectivity = np.vstack((elem_connec1, elem_connec2)).T

  # read the mesh data structure
  mesh = Mesh1D(nodes, connectivity)
  # mesh.refine_mesh(1)
  elements2node = mesh.get_mesh_coordinates(
  )    # dict: elements number with nodes coodinates

  # ====================== Pysical Problem ======================
  # define the materials
  # given JCA porous material properties
  phi = 0.99    # porosity
  sigma = 1.0567e4    # resistivity
  alpha = 1.2    # Tortuosity
  Lambda_prime = 490e-6    # Viscous characteristic length
  Lambda = 240e-6    #
  rho_1 = 9.2
  nu = 0.285
  E = 3.155e5
  eta = 0.032
  xfm = PoroElasticMaterial('xfm', phi, sigma, alpha, Lambda_prime, Lambda,
                            rho_1, E, nu, eta)

  air = Air('classical air')
  # Harmonic Acoustic problem define the frequency
  freq = 2000
  omega = 2 * np.pi * freq    # angular frequency

  k_0 = omega / Air.c
  theta = 0    #indidence angle
  ky = k_0 * np.sin(theta * np.pi / 180)
  kx = k_0 * np.cos(theta * np.pi / 180)

  # define the subdomains: domain name (material) and the elements in the domain
  air_elements = np.arange(0, int(num_elem / 2))
  xfm_elements = np.arange(int(num_elem / 2), num_elem)
  mesh.set_subdomains({air: air_elements, xfm: xfm_elements})
  check_material_compability(mesh.subdomains)
  # print(elements_set)

  order = 3    # global order of the bases
  # applied the basis on each element
  for mat, elems in mesh.subdomains.items():
    if mat.TYPE == 'Fluid':
      Pf_bases = [
          Lobbato1DElement('Pf', order, elements2node[elem]) for elem in elems
      ]    # basis for pressure in fluid domain
    elif mat.TYPE == 'Poroelastic':
      Pb_bases = [
          Lobbato1DElement('Pb', order, elements2node[elem]) for elem in elems
      ]    # basis for pressure in porous domain
      Ux_bases = [
          Lobbato1DElement('Ux', order, elements2node[elem]) for elem in elems
      ]    # basis for solid displacement in porous domain
    else:
      raise ValueError("Material type is not defined!")

  Helmholtz_fe_space = FESpace(mesh, Pf_bases)
  Biot_fe_space = FESpace(mesh, Pb_bases, Ux_bases)

  # initialize the assembler
  Helmholtz_assember = HelmholtzAssembler(Helmholtz_fe_space,
                                          mesh.subdomains,
                                          dtype=np.complex128)
  Helmholtz_assember.assembly_global_matrix(Pf_bases, 'Pf', omega)

  Biot_assember = BiotAssembler(Biot_fe_space,
                                mesh.subdomains,
                                dtype=np.complex128)
  Biot_assember.assembly_global_matrix([Pb_bases, Ux_bases], ['Pb', 'Ux'],
                                       omega)
  Assembler = CouplingAssember(mesh,
                               mesh.subdomains,
                               [Helmholtz_assember, Biot_assember],
                               coupling_type="PP_continue")
  left_hand_matrix = Assembler.assembly_gloabl_matrix()

  # import pdb; pdb.set_trace()

  # ============================= Boundary conditions =====================================
  fe_space = FESpace(mesh, Pf_bases, Pb_bases, Ux_bases)
  right_hand_vec = np.zeros(Assembler.nb_global_dofs, dtype=np.complex128)
  # essential_bcs = {'type': 'solid_displacement', 'value': P_a(-1,0), 'position': -1.}  # position is the x coordinate
  BCs_applier = ApplyBoundaryConditions(mesh, fe_space, left_hand_matrix,
                                        right_hand_vec, omega)
  # BCs_applier.apply_essential_bc(essential_bcs, var='Pf', bctype='strong')
  #  natural boundary condition
  nature_bcs_1 = {
      'type': 'total_displacement',
      'value': -u_a(-1., 0),
      'position': -1.
  }    # position is the x coordinate
  # nature_bcs_1 = {'type': 'total_displacement', 'value': 1, 'position': -1.}  # position is the x coordinate
  nature_bcs_2 = {
      'type': 'total_displacement',
      'value': u_t(1, 0),
      'position': 1.
  }    # position is the x coordinate
  nature_bcs_3 = {
      'type': 'solid_stress',
      'value': sigma_xx(1, 0),
      'position': 1.
  }    # position is the x coordinate
  BCs_applier.apply_nature_bc(nature_bcs_1, var='Pf')
  BCs_applier.apply_nature_bc(nature_bcs_2, var='Pb')
  BCs_applier.apply_nature_bc(nature_bcs_3, var='Ux')
  # import pdb; pdb.set_trace()
  # ============================= Solve the linear system ================================
  # solver the linear system
  linear_solver = LinearSolver(coupling_assember=Assembler)
  # # print("condition number:", linear_solver.condition_number(left_hand_matrix))
  # # plot_matrix_partten(left_hand_matrix)
  # import pdb;pdb.set_trace()

  linear_solver.solve(BCs_applier.left_hand_side, BCs_applier.right_hand_side)
  sol = linear_solver.u

  # ================================== Analytical Solution ======================
  # analytical solution
  xfm.set_frequency(omega)
  ana_sol_p = Fluid_Biot_Pressure(nodes, [0])
  ana_sol_u = Fluid_Biot_Displacement(nodes[int(num_elem / 2):], [0])
  # plot the solution
  post_processer_p = PostProcessField(mesh.nodes,
                                      r'1D Biot (2000$Hz$) Pressure')
  post_processer_p.plot_sol(
      (np.real(sol[:num_elem + 1]), f'FEM ($p={order}$)', 'solid'),
      (np.real(ana_sol_p[:]), 'Analytical', 'dashed'))
  # post_processer_p.plot_sol((np.real(ana_sol[:]), 'Analytical', 'solid'))
  # post_processer_p.plot_sol((np.real(sol[:num_elem+1]), 'Analytical', 'solid'))
  plt.show()
  # import pdb;pdb.set_trace()

  post_processer_u = PostProcessField(
      mesh.nodes[int(num_elem / 2):], r'1D Biot (2000$Hz$) Solid displacement')
  post_processer_u.plot_sol(
      (np.real(sol[num_elem + 1:]), f'FEM ($p=3$)', 'solid'),
      (np.real(ana_sol_u[:]), 'Analytical', 'dashed'))
  plt.show()

  error_p = post_processer_p.compute_error(sol[:num_elem + 1], ana_sol_p[:])
  error_u = post_processer_u.compute_error(sol[num_elem + 1:], ana_sol_u[:],
                                           -1)

  print("error_p: ", error_p)
  print("error_u: ", error_u)
  # if error_p < 1e-4 and error_u < 1e-3:
  #     print("Test passed!")
  #     return True
  # else:
  #     print("Test failed!")
  #     return False


if __name__ == "__main__":
  result = test_case()
