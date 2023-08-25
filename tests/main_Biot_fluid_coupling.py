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

from fem.basis import Lobbato1DElement
from fem.mesh import Mesh1D
from fem.dofhandler import DofHandler1D, GeneralDofHandler1D
from fem.physic_assembler import HelmholtzAssembler, BiotAssembler
from fem.materials import Air, Fluid, EquivalentFluid, PoroElasticMaterial
from fem.utilities import check_material_compability, display_matrix_in_array, plot_matrix_partten
from fem.solver import LinearSolver
from fem.postprocess import PostProcessField
from analytical.Biot_sol import solve_PW

def test_case():
    num_elem = 10  # number of elements
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
    elements2node = mesh.get_mesh()  # dict: elements number with nodes coodinates

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

    air = Air('classical air')
    # Harmonic Acoustic problem define the frequency
    freq = 2000
    omega = 2 * np.pi * freq  # angular frequency

    k_0 = omega/Air.c
    theta = 0  #indidence angle
    ky = k_0*np.sin(theta*np.pi/180)
    kx = k_0*np.cos(theta*np.pi/180)

    # define the subdomains: domain name (material) and the elements in the domain
    air_elements = np.arange(0, 5)
    xfm_elements = np.arange(5, num_elem)
    subdomains = {air: air_elements, xfm: xfm_elements}
    check_material_compability(subdomains)
    # print(elements_set)

    order = 3  # global order of the bases
    # applied the basis on each element
    for mat, elems in subdomains.items():
        if mat.TYPE == 'Fluid':
            Pf_bases = [Lobbato1DElement('Pf', order, elements2node[elem]) for elem in elems]  # basis for pressure in fluid domain
        elif mat.TYPE == 'Poroelastic':
            Pb_bases = [Lobbato1DElement('Pb', order, elements2node[elem]) for elem in elems]   # basis for pressure in porous domain
            Ux_bases = [Lobbato1DElement('Ux', order, elements2node[elem]) for elem in elems]  # basis for solid displacement in porous domain
        else:
            raise ValueError("Material type is not defined!")

    Helmholtz_dof_handler = GeneralDofHandler1D(('Pf'), Pf_bases)
    Biot_dof_handler = GeneralDofHandler1D(('Pb','Ux'), Pb_bases, Ux_bases)
    print(Helmholtz_dof_handler.get_num_dofs())
    print(Helmholtz_dof_handler.get_global_dofs())
    print(Biot_dof_handler.get_global_dofs())


    import time

    # initialize the assembler
    Helmholtz_assember = HelmholtzAssembler(Helmholtz_dof_handler, subdomains, dtype=np.complex128)
    Helmholtz_assember.assemble_material_K(Pf_bases, 'Pf', omega)
    Helmholtz_assember.assemble_material_M(Pf_bases, 'Pf', omega)

    K_f, M_f = Helmholtz_assember.K, Helmholtz_assember.M

    Biot_assember = BiotAssembler(Biot_dof_handler, subdomains, dtype=np.complex128)
    Biot_assember.assemble_material_K(Pb_bases, 'Pb', omega)
    Biot_assember.assemble_material_M(Pb_bases, 'Pb', omega)
    Biot_assember.assemble_material_K(Ux_bases, 'Ux', omega)
    Biot_assember.assemble_material_M(Ux_bases, 'Ux', omega)
    Biot_assember.assemble_material_C(Pb_bases, 'Ux', 'Pb', omega)

    K_p, M_p, K_u, M_u, C_pu, C_up = Biot_assember.K, Biot_assember.M, Biot_assember.K, Biot_assember.M, Biot_assember.C, Biot_assember.C.T
    
    left_hand_matrix = K_p+K_u-M_u-M_p - C_pu-C_up
    import pdb; pdb.set_trace()



    # ============================= Boundary conditions =====================================
    essential_bcs = {'type': 'solid_displacement', 'value': 0, 'position': 0.}  # position is the x coordinate
    left_hand_matrix =assembler.apply_essential_bc(left_hand_matrix, essential_bcs, var='Ux', bctype='nitsche')
    #  natural boundary condition   
    nature_bcs = {'type': 'total_displacement', 'value': 1, 'position': -1.}  # position is the x coordinate
    right_hand_vector = assembler.apply_nature_bc(nature_bcs, var='P')

    # ============================= Solve the linear system ================================
    # solver the linear system
    linear_solver = LinearSolver(dof_handler)
    # print("condition number:", linear_solver.condition_number(left_hand_matrix))
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