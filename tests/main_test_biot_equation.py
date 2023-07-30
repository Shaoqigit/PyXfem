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
from fem.dofhandler import DofHandler1D, DofHandler1DMutipleVariable
from fem.assembly import Assembler, Assembler4Biot
from fem.materials import Air, Fluid, EquivalentFluid, PoroElasticMaterial
from fem.utilities import check_material_compability, display_matrix_in_array, plot_matrix_partten
from fem.solver import LinearSolver
from fem.postprocess import PostProcessField
from fem.analytical_sol import DoubleleLayerKundltTube


def test_case():
    num_elem = 100  # number of elements
    num_nodes = num_elem + 1  # number of nodes

    nodes = np.linspace(0, 1, num_nodes)

    elem_connec1 = np.arange(0, num_elem)
    elem_connec2 = np.arange(1, num_nodes)
    connectivity = np.vstack((elem_connec1, elem_connec2)).T
    # print(connectivity)


    # read the mesh data structure
    mesh = Mesh1D(nodes, connectivity)
    # mesh.refine_mesh(1)

    elements_set = mesh.get_mesh()  # dict: elements number with nodes coodinates
    # print(elements_set)

    P_bases = []  # basis applied on each element, could be different order and type
    Ux_bases = []
    order = 2  # global order of the bases
    # applied the basis on each element
    for key, elem in elements_set.items():
        P_basis = Lobbato1DElement('P', order, elem)
        Ux_basis = Lobbato1DElement('Ux', order, elem)
        P_bases.append(P_basis)
        Ux_bases.append(Ux_basis)
        # print(basis.ke)
        # print(basis.me)

    # handler the dofs: map the basis to mesh
    dof_handler = DofHandler1DMutipleVariable(mesh, P_bases, Ux_bases)
    print(dof_handler.get_num_dofs())

    print("global dofs index: ", dof_handler.get_global_dofs())
    print("bases global dofs index: ", dof_handler.base4global_dofs())
    print("global dofs index for P: ", dof_handler.get_global_dofs_by_base('P'))
    print("global dofs index for Ux: ", dof_handler.get_global_dofs_by_base('Ux'))

    # print(dof_handler.num_internal_dofs)
    # print(dof_handler.num_external_dofs)

    # ====================== Pysical Problem ======================
    # define the materials
    # given JCA porous material properties
    phi          = 0.98  # porosity
    sigma        = 3.75e3 # resistivity
    alpha        = 1.17  # Tortuosity
    Lambda_prime = 742e-6  # Viscous characteristic length
    Lambda      = 110e-6  # 
    rho_1 = 22.1
    nu = 0.39
    E = 70e3
    eta = 0.265
    xfm = PoroElasticMaterial('xfm', phi, sigma, alpha, Lambda_prime, Lambda, rho_1, E, nu, eta)

    # Harmonic Acoustic problem define the frequency
    freq = 2000
    omega = 2 * np.pi * freq  # angular frequency

    # define the subdomains: domain name (material) and the elements in the domain
    xfm_elements = np.arange(0, num_nodes)
    subdomains = {xfm: xfm_elements}
    check_material_compability(subdomains)


    # initialize the assembler
    assembler = Assembler4Biot(dof_handler, subdomains, dtype=np.complex128)

    K_p= assembler.assemble_material_K(P_bases, 'P', omega)  # global stiffness matrix with material attribution
    M_p= assembler.assemble_material_M(P_bases, 'P', omega)  # global mass matrix with material attribution
    # print("K_g:", assembler.get_matrix_in_array(K_g))
    # print("M_g:", assembler.get_matrix_in_array(M_g))

    K_u= assembler.assemble_material_K(Ux_bases, 'Ux', omega)  # global stiffness matrix with material attribution
    M_u= assembler.assemble_material_M(Ux_bases, 'Ux', omega)  # global mass matrix with material attribution

    C_pu= assembler.assemble_material_C(P_bases, 'P', 'Ux', omega)  # global coupling matrix with material attribution
    C_up= C_pu.T  # global coupling matrix with material attribution 
    # construct linear system
    # K_u-M_u - C+(K_p-M_p -C.T)
    left_hand_matrix = K_p+K_u-M_u-M_p - C_pu-C_up
    plot_matrix_partten(left_hand_matrix)

    essential_bcs = {'type': 'solid_displacement', 'value': 0, 'position': num_elem}
    assembler.apply_essential_bc(left_hand_matrix, essential_bcs, var='Ux',)
    # print(assembler.get_matrix_in_array(left_hand_matrix))
    #  natural boundary condition   
    nature_bcs = {'type': 'total_displacement', 'value': 1*np.exp(-1j*omega), 'position': 0}
    right_hand_vector = assembler.apply_nature_bc(nature_bcs)
    # print(right_hand_vector)

    # solver the linear system
    linear_solver = LinearSolver(dof_handler)
    # plot_matrix_partten(left_hand_matrix)
    linear_solver.solve(left_hand_matrix, right_hand_vector)
    sol = linear_solver.u


    # ====================== Analytical Solution ======================
    # analytical solution

    import pdb; pdb.set_trace()
    # plot the solution
    post_processer = PostProcessField(mesh.nodes, r'1D Helmholtz (2000$Hz$)')
    # post_processer.plot_sol((np.real(sol), f'FEM ($p=3$)', 'solid'), (np.real(ana_sol), 'Analytical', 'dashed'))
    post_processer.plot_sol((np.real(sol[:101]), f'FEM ($p=3$)', 'solid'))
    plt.show()

    # error = post_processer.compute_error(sol, ana_sol)
    # print("error:", error)
    # if error < 1e-4:
    #     print("Test passed!")
    #     return True
    # else:
    #     print("Test failed!")
    #     return False

if __name__ == "__main__":
    result = test_case()