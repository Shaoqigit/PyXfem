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

from SAcouS.acxfem.basis import Lobbato1DElement
from SAcouS.acxfem.mesh import Mesh1D
from SAcouS.acxfem.dofhandler import FESpace
from SAcouS.acxfem.assembly import Assembler
from SAcouS.acxfem.materials import Air, EquivalentFluid
from SAcouS.acxfem.physic_assembler import HelmholtzAssembler
from SAcouS.acxfem.utilities import check_material_compability
from SAcouS.acxfem.solver import LinearSolver

from SAcouS.acxmor.modal_reduction import EigenSolver, ModalReduction
from SAcouS.acxfem.postprocess import PostProcessField, PostProcessFRF
from analytical.fluid_sol import ImpedenceKundltTube


def FEM_model(omega, assembler, nature_bcs):
    assembler.initial_matrix()
    K_g= assembler.assemble_material_K(omega)  # global stiffness matrix with material attribution
    M_g= assembler.assemble_material_M(omega)  # global mass matrix with material attribution
    right_hand_vector = assembler.assemble_nature_bc(nature_bcs)
    left_hand_matrix = K_g-M_g

    # solver the linear system
    linear_solver = LinearSolver(fe_space)
    # import pdb; pdb.set_trace() 
    # plot_matrix_partten(left_hand_matrix)
    linear_solver.solve(left_hand_matrix, right_hand_vector)
    sol = linear_solver.u

    return sol


def reduction_model(omega, K_r, M_r, f_r):

    left_hand_matrix = K_r-omega**2*M_r
    # solver the linear system
    reduced_sol = modal_reduction_method.solve(left_hand_matrix, f_r)
    
    sol = modal_reduction_method.recover_sol(reduced_sol)

    return sol
    

if __name__ == "__main__":
    # ====================== Pysical Problem ======================
    # define the materials
    air = Air('classical air')

    # given JCA porous material properties
    phi          = 0.98  # porosity
    sigma        = 3.75e3 # resistivity
    alpha        = 1.17  # Tortuosity
    Lambda_prime = 742e-6  # Viscous characteristic length
    Lambda      = 110e-6  # 
    xfm = EquivalentFluid('xfm', phi, sigma, alpha, Lambda_prime, Lambda)

    # Harmonic Acoustic problem define the frequency
    freq = np.linspace(10, 2000, 100)  # frequency range
    omegas = 2 * np.pi * freq  # angular frequency

    num_elem = 500  # number of elements
    num_nodes = num_elem + 1  # number of nodes
    nodes = np.linspace(-1, 0, num_nodes)
    elem_connec1 = np.arange(0, num_elem)
    elem_connec2 = np.arange(1, num_nodes)
    connectivity = np.vstack((elem_connec1, elem_connec2)).T
    # read the mesh data structure
    mesh = Mesh1D(nodes, connectivity)
    elements2node = mesh.get_mesh()  # dict: elements number with nodes coodinates

    # define the subdomains: domain name (material) and the elements in the domain
    air_elements = np.arange(0, int(num_elem/2))
    xfm_elements = np.arange(int(num_elem/2), num_elem)
    subdomains = {air: air_elements, xfm: xfm_elements}
    check_material_compability(subdomains)

    order = 2  # global order of the bases
    # applied the basis on each element
    Pf_bases = []
    for mat, elems in subdomains.items():
        if mat.TYPE == 'Fluid':
            Pf_bases += [Lobbato1DElement('Pf', order, elements2node[elem]) for elem in elems] 
        # print(basis.ke)

    # handler the dofs: map the basis to mesh
    fe_space = FESpace(mesh, subdomains, Pf_bases)

    # initialize the assembler
    Helmholtz_assember = HelmholtzAssembler(fe_space, subdomains, dtype=np.complex128)
    K_w = Helmholtz_assember.assemble_material_K(Pf_bases, 'Pf', omega=1)
    M_w = Helmholtz_assember.assemble_material_M(Pf_bases, 'Pf', omega=1)

    nature_bcs = {'type': 'total_displacement', 'value': 1, 'position': 0}
    impedence_bcs = {'type': 'impedence', 'value': 0., 'position': num_elem}

    import time
    start = time.time()
    sol_fem = np.zeros(100, dtype=np.complex128)
    for i, omega in enumerate(omegas):
        sol_fem[i] = FEM_model(omega, Helmholtz_assember, nature_bcs)[int(num_nodes/2)]
    end = time.time()
    print("FEM solving time: ", end-start)

    # ====================== Modal Reduction ======================
    start = time.time()
    nb_modes = 50
    Helmholtz_assember.initial_matrix()
    # import pdb; pdb.set_trace()
    K_w= Helmholtz_assember.assemble_material_K(omega=1)  # global stiffness matrix no frequency dependent material
    M_w= Helmholtz_assember.assemble_material_M(omega=1)  # global mass matrix: the eigen freq should be omega^2
    # construct modal domain
    eigen_value_solver = EigenSolver(fe_space)
    eig_omega_sq, modes = eigen_value_solver.solve(K_w, M_w, nb_modes)
    
    modal_reduction_method = ModalReduction(K_w, M_w, modes)

    #  natural boundary condition   
    right_hand_vector = Helmholtz_assember.assemble_nature_bc(nature_bcs)
    # ====================== Reduced System ======================
    # reduced the system
    K_r = modal_reduction_method.projection(K_w)
    M_r = modal_reduction_method.projection(M_w)
    f_r = modal_reduction_method.projection(right_hand_vector)

    sol_modal = np.zeros(100, dtype=np.complex128)
    for i, omega in enumerate(omegas):
        sol_modal[i] = reduction_model(omega, K_r, M_r, f_r)[int(num_nodes/2)]
    end = time.time()
    print("Modal reduction solving time: ", end-start)
        # ====================== Analytical Solution ======================
    # analytical solution
    # kundlt_tube = ImpedenceKundltTube(mesh, air, omega, nature_bcs, impedence_bcs)
    # ana_sol = np.zeros(num_nodes, dtype=np.complex128)  #initialize the analytical solution vector
    # kundlt_tube.sol_on_nodes(ana_sol, sol_type='pressure')

    # # plot the solution
    post_process = PostProcessFRF(freq, r'1D Helmholtz FRF', 'SPL(dB)')
    post_process.plot_sol((np.real(sol_fem), f'FEM (dofs$={num_elem}$)', 'solid'), (np.real(sol_modal), f'Modal reduction ($m={nb_modes}$)', 'dashed'))
    plt.show(block=False)
    # plt.show()
    plt.pause(1)
    plt.close('all')

    # compute the error
    error = post_process.compute_error(sol_fem, sol_modal)
    print("error: ", error)
    if error < 1e-4:
        print("Test passed!")
