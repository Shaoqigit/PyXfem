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

# first tuturial, in which almost all the modules are used to demonstrate the usage of the PyAcoustiX
# the problem is a 1D Helmholtz problem with two materials: air and JCA porous material
# the analytical solution is given by Kundlt tube as the reference solution to compare with the FEM solution

import sys
sys.path.append('/home/shaoqi/Devlop/PyXfem/PyAcoustiX/')
import numpy as np
import matplotlib.pyplot as plt

from SAcouS.acxfem.basis import Lobbato1DElement
from SAcouS.acxfem.mesh import Mesh1D
from SAcouS.acxfem.dofhandler import GeneralDofHandler1D, FESpace
from SAcouS.acxfem.physic_assembler import HelmholtzAssembler
from SAcouS.acxfem.BCs_impose import ApplyBoundaryConditions
from SAcouS.acxfem.materials import Air, EquivalentFluid
from SAcouS.acxfem.utilities import check_material_compability, plot_matrix_partten
from SAcouS.acxfem.solver import LinearSolver
from SAcouS.acxfem.postprocess import PostProcessField
from analytical.fluid_sol import DoubleleLayerKundltTube


def tutor_case_1():
    # ====================== Pysical Problem ======================
    print("=====================set up the physical problem===========================")
    # define the materials
    # classical air as propagation medium
    air = Air('classical air')
    # given JCA porous material properties
    phi          = 0.98  # porosity
    sigma        = 3.75e3 # resistivity
    alpha        = 1.17  # Tortuosity
    Lambda_prime = 742e-6  # Viscous characteristic length
    Lambda      = 110e-6  # 
    xfm = EquivalentFluid('xfm', phi, sigma, alpha, Lambda_prime, Lambda)  # equivalent fluid model for porous material

    # Harmonic Acoustic problem define the frequency
    freq = 2000
    omega = 2 * np.pi * freq  # angular frequency

    # ====================== Mesh and basis definition ======================
    print("=====================set up the mesh===========================")
    num_elem = 100  # number of elements
    num_nodes = num_elem + 1  # number of nodes
    nodes = np.linspace(-1, 1, num_nodes)

    elem_connec1 = np.arange(0, num_elem)
    elem_connec2 = np.arange(1, num_nodes)
    connectivity = np.vstack((elem_connec1, elem_connec2)).T
    # crate 1D mesh with nodes coordinates and connectivity
    mesh = Mesh1D(nodes, connectivity)
    elements2node = mesh.get_mesh()  # dict: elements number with nodes coodinates
    print(f"The of {num_elem} elements mesh is created.")

    # define the subdomains: domain name (material) and the elements in the domain
    air_elements = np.arange(0, int(num_elem/2))
    xfm_elements = np.arange(int(num_elem/2), num_elem)
    subdomains = {air: air_elements, xfm: xfm_elements}

    # check if the materials defined in the mesh domains are compatible
    check_material_compability(subdomains)

    print("=====================set up the FEM bases===========================")
    order = 2  # global order of the bases, which can be defined respectively on each elements
    # applied the basis on each element
    Pf_bases = []
    for mat, elems in subdomains.items():
        if mat.TYPE == 'Fluid':
            Pf_bases += [Lobbato1DElement('Pf', order, elements2node[elem]) for elem in elems]  # arguments: name, order, nodes coordinates

    # handler the dofs: arange the dofs numbering for current bases
    Helmholtz_dof_handler = GeneralDofHandler1D(['Pf'], Pf_bases)  # args: bases name, bases
    print(f"The number of global dofs is {Helmholtz_dof_handler.get_nb_dofs()}.")
    # initialize the assembler
    # here two fluid-like materials are defined, so the assembler will be initialized with high wraped HelmholtzAssembler (level-wrapped general assembly can be used as well, but more complicated)
    Helmholtz_assember = HelmholtzAssembler(Helmholtz_dof_handler, subdomains, dtype=np.complex128)  # args: dof_handler, subdomains, dtype (type of left hand matrix and right hand vector)
    # assembly the global matrix with omega
    Helmholtz_assember.assembly_global_matrix(Pf_bases, 'Pf', omega) 
    left_hand_matrix = Helmholtz_assember.get_global_matrix()  # get the global matrix
    plot_matrix_partten(left_hand_matrix)  # display the matrix sparse partten


    # ====================== Boundary Conditions ======================
    print("=====================set up the boundary conditions===========================")
    fe_space = FESpace(mesh, subdomains, Pf_bases)
    right_hand_vec = np.zeros(Helmholtz_assember.nb_global_dofs, dtype=np.complex128)

    #  natural boundary condition   
    nature_bcs = {'type': 'velocity', 'value': 1*np.exp(-1j*omega), 'position': -1}
    BCs_applier = ApplyBoundaryConditions(mesh, fe_space, left_hand_matrix, right_hand_vec, omega)
    BCs_applier.apply_nature_bc(nature_bcs, var='Pf')
    print(f"{nature_bcs['type']} boundary condition is applied at {nature_bcs['position']}.")

    print("=====================solve the linear system===========================")
    # solver the linear system
    linear_solver = LinearSolver(dof_handler=Helmholtz_dof_handler)
    linear_solver.solve(left_hand_matrix, right_hand_vec)
    sol = linear_solver.u


    # ====================== Analytical Solution ======================
    # analytical solution
    kundlt_tube = DoubleleLayerKundltTube(mesh, air, xfm, omega, nature_bcs)
    ana_sol = np.zeros(num_nodes, dtype=np.complex128)  #initialize the analytical solution vector
    kundlt_tube.sol_on_nodes(ana_sol, sol_type='pressure')

    # plot the solution
    post_processer = PostProcessField(mesh.nodes, r'1D Helmholtz (2000$Hz$)')
    post_processer.plot_sol((np.real(sol), f'FEM ($p=3$)', 'solid'), (np.real(ana_sol), 'Analytical', 'dashed'))
    plt.show()


if __name__ == "__main__":
    print("==================== AcoustiX 1D is running===========================")
    result = tutor_case_1()