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
import os
current_dir = os.path.dirname(os.path.realpath(__file__))
working_dir = os.path.join(current_dir , "..")
import sys
sys.path.append(working_dir)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import spy

from fem.basis import Lobbato1DElement
from fem.mesh import Mesh1D
from fem.dofhandler import DofHandler1D
from fem.assembly import Assembler
from fem.materials import Air, Fluid, EquivalentFluid
from fem.utilities import check_material_compability, display_matrix_in_array, plot_matrix_partten
from fem.solver import LinearSolver

from mor.modal_reduction import EigenSolver, ModalReduction
from fem.postprocess import PostProcessField
from analytical.fluid_sol import ImpedenceKundltTube


def FEM_model(omega, assembler, nature_bcs):
    K_g= assembler.assemble_material_K(omega)  # global stiffness matrix with material attribution
    M_g= assembler.assemble_material_M(omega)  # global mass matrix with material attribution

    right_hand_vector = assembler.assemble_nature_bc(nature_bcs)

    left_hand_matrix = K_g-M_g

    # solver the linear system
    linear_solver = LinearSolver(dof_handler)
    # plot_matrix_partten(left_hand_matrix)
    linear_solver.solve(left_hand_matrix, right_hand_vector)
    sol = linear_solver.u

    return sol


def reduction_model(omega, assembler, nature_bcs, nb_modes):

    K_w= assembler.assemble_material_K(omega=1)  # global stiffness matrix no frequency dependent material
    M_w= assembler.assemble_material_M(omega=1)  # global mass matrix: the eigen freq should be omega^2
  
    # construct modal domain
    import pdb; pdb.set_trace()
    eigen_value_solver = EigenSolver(dof_handler)
    eig_omega_sq, modes = eigen_value_solver.solve(K_w, M_w, nb_modes)
    eig_freqs = np.sqrt(eig_omega_sq)/(2*np.pi)
    
    modal_reduction_method = ModalReduction(eig_freqs, modes)

    # print(assembler.get_matrix_in_array(left_hand_matrix))
    #  natural boundary condition   
    assembler.omega = omega
    nature_bcs = {'type': 'velocity', 'value': 1, 'position': 0}
    right_hand_vector = assembler.assemble_nature_bc(nature_bcs)
    # print(display_matrix_in_array(C_g))

    # ====================== Reduced System ======================
    # reduced the system
    K_r = modal_reduction_method.projection(K_w)
    M_r = modal_reduction_method.projection(M_w)
    f_r = modal_reduction_method.projection(right_hand_vector)

    left_hand_matrix = K_r-omega**2*M_r


    # solver the linear system
    reduced_sol = modal_reduction_method.solve(left_hand_matrix, f_r)
    
    import pdb; pdb.set_trace()
    sol = modal_reduction_method.recover_sol(reduced_sol)

    return sol
    



if __name__ == "__main__":
    # ====================== Pysical Problem ======================
    # define the materials
    air = Air('classical air')
    # Harmonic Acoustic problem define the frequency
    freq = 2000
    omega = 2 * np.pi * freq  # angular frequency

    num_elem = 400  # number of elements
    num_nodes = num_elem + 1  # number of nodes
    nodes = np.linspace(-1, 0, num_nodes)
    elem_connec1 = np.arange(0, num_elem)
    elem_connec2 = np.arange(1, num_nodes)
    connectivity = np.vstack((elem_connec1, elem_connec2)).T
    # read the mesh data structure
    mesh = Mesh1D(nodes, connectivity)
    elements_set = mesh.get_mesh()  # dict: elements number with nodes coodinates

    # define the subdomains: domain name (material) and the elements in the domain
    air_elements = np.arange(0, num_nodes)
    subdomains = {air: air_elements}
    check_material_compability(subdomains)

    bases = []  # basis applied on each element, could be different order and type
    order = 1  # global order of the bases
    # applied the basis on each element
    for key, elem in elements_set.items():
        basis = Lobbato1DElement('P', order, elem)
        bases.append(basis)

    # handler the dofs: map the basis to mesh
    dof_handler = DofHandler1D(mesh, bases)

    assembler = Assembler(dof_handler, bases, subdomains, dtype=np.complex128)
    nature_bcs = {'type': 'velocity', 'value': 1*np.exp(1j), 'position': 0}
    impedence_bcs = {'type': 'impedence', 'value': 0., 'position': num_elem}


    sol_fem = FEM_model(omega, assembler, nature_bcs)
    nb_modes = 100
    sol_modal = reduction_model(omega, assembler, nature_bcs, nb_modes)

        # ====================== Analytical Solution ======================
    # analytical solution
    kundlt_tube = ImpedenceKundltTube(mesh, air, omega, nature_bcs, impedence_bcs)
    ana_sol = np.zeros(num_nodes, dtype=np.complex128)  #initialize the analytical solution vector
    kundlt_tube.sol_on_nodes(ana_sol, sol_type='pressure')

    # # plot the solution
    post_process = PostProcessField(mesh.nodes, r'1D Helmholtz (2000$Hz$)')
    post_process.plot_sol((np.real(sol_fem), f'FEM ($p=1$)', 'solid'), (np.real(sol_modal), f'Modal reduction ($m={nb_modes}$)', 'solid'), (np.real(ana_sol), 'Analytical', 'dashed'))
    plt.show()
    # plt.pause(1)
    # plt.close('all')

    # compute the error
    error = post_process.compute_error(sol_fem, sol_modal)
    print("error: ", error)
