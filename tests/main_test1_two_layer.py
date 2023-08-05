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
import meshio
import matplotlib.pyplot as plt
from matplotlib.pyplot import spy

from fem.basis import Lobbato1DElement
from fem.mesh import Mesh1D
from fem.dofhandler import DofHandler1D
from fem.assembly import Assembler
from fem.materials import Air, Fluid, EquivalentFluid
from fem.utilities import check_material_compability, display_matrix_in_array, plot_matrix_partten
from fem.solver import LinearSolver
from fem.postprocess import PostProcessField
from analytical.fluid_sol import DoubleleLayerKundltTube


def test_case_1():
    num_elem = 1000  # number of elements
    num_nodes = num_elem + 1  # number of nodes

    nodes = np.linspace(-1, 1, num_nodes)

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
    order = 1  # global order of the bases
    # applied the basis on each element
    for key, elem in elements_set.items():
        basis = Lobbato1DElement('P', order, elem)
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
    water=Fluid('water', 997, 1481)

    # given JCA porous material properties
    phi          = 0.98  # porosity
    sigma        = 3.75e3 # resistivity
    alpha        = 1.17  # Tortuosity
    Lambda_prime = 742e-6  # Viscous characteristic length
    Lambda      = 110e-6  # 
    xfm = EquivalentFluid('xfm', phi, sigma, alpha, Lambda_prime, Lambda)

    # Harmonic Acoustic problem define the frequency
    freq = 8000
    omega = 2 * np.pi * freq  # angular frequency

    # define the subdomains: domain name (material) and the elements in the domain
    air_elements = np.arange(0, num_elem/2)
    xfm_elements = np.arange(num_elem/2, num_nodes)
    subdomains = {air: air_elements, xfm: xfm_elements}
    check_material_compability(subdomains)


    # initialize the assembler
    assembler = Assembler(dof_handler, bases, subdomains, dtype=np.complex128)

    K_g= assembler.assemble_material_K(omega)  # global stiffness matrix with material attribution
    M_g= assembler.assemble_material_M(omega)  # global mass matrix with material attribution
    # print("K_g:", assembler.get_matrix_in_array(K_g))
    # print("M_g:", assembler.get_matrix_in_array(M_g))

    # construct linear system
    left_hand_matrix = K_g-M_g
    # plot_matrix_partten(left_hand_matrix)

    # print(assembler.get_matrix_in_array(left_hand_matrix))
    #  natural boundary condition   
    nature_bcs = {'type': 'velocity', 'value': 1*np.exp(-1j*omega), 'position': 0}
    right_hand_vector = assembler.assemble_nature_bc(nature_bcs)
    # print(right_hand_vector)

    # solver the linear system
    linear_solver = LinearSolver(dof_handler)
    # left_hand_matrix, right_hand_vector = linear_solver.optimize_matrix_pattern(left_hand_matrix, right_hand_vector)
    # left_hand_matrix = left_hand_matrix[optimal_order, optimal_order]
    # print(left_hand_matrix[0], left_hand_matrix[1])
    # plot_matrix_partten(left_hand_matrix)
    linear_solver.solve(left_hand_matrix, right_hand_vector)
    sol = linear_solver.u
    # print("solution:", abs(sol))


    # ====================== Analytical Solution ======================
    # analytical solution
    kundlt_tube = DoubleleLayerKundltTube(mesh, air, xfm, omega, nature_bcs)
    ana_sol = np.zeros(num_nodes, dtype=np.complex128)  #initialize the analytical solution vector
    kundlt_tube.sol_on_nodes(ana_sol, sol_type='pressure')

    # plot the solution
    post_processer = PostProcessField(mesh.nodes, r'1D Helmholtz (2000$Hz$)')
    post_processer.plot_sol((np.real(sol), f'FEM ($p=3$)', 'solid'), (np.real(ana_sol), 'Analytical', 'dashed'))
    plt.show()

    error = post_processer.compute_error(sol, ana_sol)
    print("error:", error)
    if error < 1e-5:
        print("Test passed!")
        return True
    else:
        print("Test failed!")
        return False

if __name__ == "__main__":
    result = test_case_1()