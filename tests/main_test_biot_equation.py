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
from analytical.Biot_sol import solve_PW

def test_case():
    num_elem = 1000  # number of elements
    num_nodes = num_elem + 1  # number of nodes

    nodes = np.linspace(-1, 0, num_nodes)
    # print()
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
        Ux_basis = Lobbato1DElement('Ux', order, elem)
        P_basis = Lobbato1DElement('P', order, elem)
        P_bases.append(P_basis)
        Ux_bases.append(Ux_basis)

    # handler the dofs: map the basis to mesh
    dof_handler = DofHandler1DMutipleVariable(mesh, Ux_bases, P_bases)
    print(dof_handler.get_num_dofs())



    # ====================== Pysical Problem ======================
    # define the materials
    # given JCA porous material properties
    phi          = 0.99  # porosity
    sigma        = 1.0567e4 # resistivity
    alpha        = 1.2  # Tortuosity
    Lambda_prime = 490e-6  # Viscous characteristic length
    Lambda      = 240e-6  # 
    rho_1 = 9.2
    nu = 0.285
    E = 3.155e5
    eta = 0.032
    xfm = PoroElasticMaterial('xfm', phi, sigma, alpha, Lambda_prime, Lambda, rho_1, E, nu, eta)

    # Harmonic Acoustic problem define the frequency
    freq = 2000
    omega = 2 * np.pi * freq  # angular frequency

    k_0 = omega/Air.c
    theta = 0  #indidence angle
    ky = k_0*np.sin(theta*np.pi/180)
    kx = k_0*np.cos(theta*np.pi/180)

    # define the subdomains: domain name (material) and the elements in the domain
    xfm_elements = np.arange(0, num_nodes)
    subdomains = {xfm: xfm_elements}
    check_material_compability(subdomains)

    import time
    # initialize the assembler
    assembler = Assembler4Biot(dof_handler, subdomains, dtype=np.complex128)

    start = time.perf_counter()
    K_p= assembler.assemble_material_K(P_bases, 'Ux', omega)  # global stiffness matrix with material attribution
    end = time.perf_counter()
    print("Elapsed (with compilation) = {}s".format((end - start)))

    
    M_p= assembler.assemble_material_M(P_bases, 'Ux', omega)  # global mass matrix with material attribution
    assembler.initial_matrix()
    K_u= assembler.assemble_material_K(Ux_bases, 'P', omega)  # global stiffness matrix with material attribution
    M_u= assembler.assemble_material_M(Ux_bases, 'P', omega)  # global mass matrix with material attribution

    C_up= assembler.assemble_material_C(P_bases, 'Ux', 'P', omega)  # global coupling matrix with material attribution
    
    C_pu= C_up.T  # global coupling matrix with material attribution 
    # construct linear system
    left_hand_matrix = K_p+K_u-M_u-M_p - C_pu-C_up

    # ============================= Boundary conditions =====================================
    essential_bcs = {'type': 'solid_displacement', 'value': 0, 'position': 0.}  # position is the x coordinate
    left_hand_matrix =assembler.apply_essential_bc(left_hand_matrix, essential_bcs, var='Ux', bctype='nitsche')
    #  natural boundary condition   
    nature_bcs = {'type': 'total_displacement', 'value': 1, 'position': -1.}  # position is the x coordinate
    right_hand_vector = assembler.apply_nature_bc(nature_bcs, var='P')

    # solver the linear system
    linear_solver = LinearSolver(dof_handler)
    # plot_matrix_partten(left_hand_matrix)
    linear_solver.solve(left_hand_matrix, right_hand_vector)
    sol = linear_solver.u


    # ================================== Analytical Solution ======================
    # analytical solution
    xfm.set_frequency(omega)
    ana_sol = solve_PW(xfm,ky,nodes,1)
    # plot the solution
    post_processer_p = PostProcessField(mesh.nodes, r'1D Biot (2000$Hz$) Pressure')
    post_processer_p.plot_sol((np.real(sol[num_elem+1:]), f'FEM ($p=3$)', 'solid'), (np.real(ana_sol[4,:]), 'Analytical', 'dashed'))
    # post_processer.plot_sol((np.real(sol[:101]), f'FEM ($p=3$)', 'solid'))
    # plt.show()

    post_processer_u = PostProcessField(mesh.nodes, r'1D Biot (2000$Hz$) Solid displacement')
    post_processer_u.plot_sol((np.real(sol[:num_elem+1]), f'FEM ($p=3$)', 'solid'), (np.real(ana_sol[1,:]), 'Analytical', 'dashed'))
    # post_processer.plot_sol((np.real(sol[:101]), f'FEM ($p=3$)', 'solid'))
    # plt.show()
    
    import pdb; pdb.set_trace()

    error_p = post_processer_p.compute_error(sol[num_elem+1:], ana_sol[4,:])
    error_u = post_processer_u.compute_error(sol[:num_elem+1], ana_sol[1,:], -1)

    print("error_p: ", error_p)
    print("error_u: ", error_u)
    if error_p < 1e-4 and error_u < 1e-3:
        print("Test passed!")
        return True
    else:
        print("Test failed!")
        return False

if __name__ == "__main__":
    result = test_case()