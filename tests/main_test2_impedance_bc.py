# This file is part of PyXfem, a software distributed under the MIT license.
# For any question, please contact the authors cited below.
#
# Copyright (c) 2023
# 	Shaoqi WU <shaoqiwu@outlook.com
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

from src.basis import Lobbato1DElement
from src.mesh import Mesh1D
from src.dofhandler import DofHandler1D
from src.assembly import Assembler
from src.materials import Air, Fluid, EquivalentFluid
from src.utilities import check_material_compability, display_matrix_in_array, plot_matrix_partten
from src.solver import LinearSolver
from src.postprocess import PostProcessField
from src.analytical_sol import ImpedenceKundltTube


def test_case_2():
    num_elem = 100  # number of elements
    num_nodes = num_elem + 1  # number of nodes

    nodes = np.linspace(-1, 0, num_nodes)

    elem_connec1 = np.arange(0, num_elem)
    elem_connec2 = np.arange(1, num_nodes)
    connectivity = np.vstack((elem_connec1, elem_connec2)).T
    # print(connectivity)


    # read the mesh data structure
    mesh = Mesh1D(nodes, connectivity)
    # mesh.refine_mesh(1)

    elements_set = mesh.get_mesh()  # dict: elements number with nodes coodinates
    # print(elements_set)

    bases = []  # basis applied on each element, could be different order and type
    order = 4  # global order of the bases
    # applied the basis on each element
    for key, elem in elements_set.items():
        basis = Lobbato1DElement(order, elem)
        bases.append(basis)
        # print(basis.ke)
        # print(basis.me)

    # handler the dofs: map the basis to mesh
    dof_handler = DofHandler1D(mesh, bases)
    # print("global dofs index: ", dof_handler.get_global_dofs())
    # print(dof_handler.get_num_dofs())
    # print(dof_handler.num_internal_dofs)
    # print(dof_handler.num_external_dofs)

    # ====================== Pysical Problem ======================
    # define the materials
    air = Air('classical air')

    # Harmonic Acoustic problem define the frequency
    freq = 2000
    omega = 2 * np.pi * freq  # angular frequency

    # define the subdomains: domain name (material) and the elements in the domain
    air_elements = np.arange(0, num_nodes)
    subdomains = {air: air_elements}
    check_material_compability(subdomains)


    # initialize the assembler
    assembler = Assembler(dof_handler, bases, subdomains, dtype=np.complex128)

    K_g= assembler.assemble_material_K(omega)  # global stiffness matrix with material attribution
    M_g= assembler.assemble_material_M(omega)  # global mass matrix with material attribution

    # construct linear system

    # print(assembler.get_matrix_in_array(left_hand_matrix))
    #  natural boundary condition   
    nature_bcs = {'type': 'velocity', 'value': 1, 'position': 0}
    impedence_bcs = {'type': 'impedence', 'value': 0.2, 'position': num_elem}
    right_hand_vector = assembler.assemble_nature_bc(nature_bcs)
    C_g = assembler.assemble_impedance_bc(impedence_bcs)
    # print(display_matrix_in_array(C_g))

    left_hand_matrix = K_g-M_g+C_g


    # solver the linear system
    linear_solver = LinearSolver(dof_handler)
    # plot_matrix_partten(left_hand_matrix)
    linear_solver.solve(left_hand_matrix, right_hand_vector)
    sol = linear_solver.u
    # print("solution:", abs(sol))


    # ====================== Analytical Solution ======================
    # analytical solution
    kundlt_tube = ImpedenceKundltTube(mesh, air, omega, nature_bcs, impedence_bcs)
    ana_sol = np.zeros(num_nodes, dtype=np.complex128)  #initialize the analytical solution vector
    kundlt_tube.sol_on_nodes(ana_sol, sol_type='pressure')

    # # plot the solution
    post_process = PostProcessField(mesh.nodes, r'1D Helmholtz (2000$Hz$)')
    post_process.plot_sol((np.real(sol), f'FEM ($p=3$)', 'solid'), (np.real(ana_sol), 'Analytical', 'dashed'))
    plt.show()

    # compute the error
    error = post_process.compute_error(sol, ana_sol)
    # print("error: ", error)
    if error < 1e-4:
        print("Test passed!")
        return True
    else:
        print("Test failed!")
        return False

if __name__ == "__main__":
    result = test_case_2()