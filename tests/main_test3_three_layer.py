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
from matplotlib import animation
from matplotlib.pyplot import spy

from SAcouS.acxfem.basis import Lobbato1DElement
from SAcouS.acxfem.mesh import Mesh1D
from SAcouS.acxfem.dofhandler import DofHandler1D
from SAcouS.acxfem.assembly import Assembler
from SAcouS.acxfem.materials import Air, Fluid, EquivalentFluid
from SAcouS.acxfem.utilities import check_material_compability, display_matrix_in_array, plot_matrix_partten
from SAcouS.acxfem.solver import LinearSolver
from SAcouS.acxfem.postprocess import PostProcessField
from analytical.fluid_sol import DoubleleLayerKundltTube



num_elem = 200  # number of elements
num_nodes = num_elem + 1  # number of nodes

nodes = np.linspace(-1, 0.5, num_nodes)

elem_connec1 = np.arange(0, num_elem)
elem_connec2 = np.arange(1, num_nodes)
connectivity = np.vstack((elem_connec1, elem_connec2)).T
# print(connectivity)


# read the mesh data structure
mesh = Mesh1D(nodes, connectivity)
# mesh.refine_mesh(1)

elements_set = mesh.get_mesh()  # dict: elements number with nodes coodinates
print(elements_set)

bases = []  # basis applied on each element, could be different order and type
order = 2  # global order of the bases
# applied the basis on each element
for key, elem in elements_set.items():
    basis = Lobbato1DElement('Pf', order, elem)
    bases.append(basis)
    # print(basis.ke)
    # print(basis.me)

# handler the dofs: map the basis to mesh
dof_handler = DofHandler1D(mesh, bases)
# print("global dofs index: ", dof_handler.get_global_dofs())

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
foam = EquivalentFluid('xfm', phi, sigma, alpha, Lambda_prime, Lambda)
# second porous layer
phi          = 0.98  # porosity
sigma        = 13.5e3 # resistivity
alpha        = 1.7  # Tortuosity
Lambda_prime = 160e-6  # Viscous characteristic length
Lambda      = 80e-6  # 
xfm = EquivalentFluid('xfm', phi, sigma, alpha, Lambda_prime, Lambda)

# Harmonic Acoustic problem define the frequency
freq = 1000
omega = 2 * np.pi * freq  # angular frequency

# define the subdomains: domain name (material) and the elements in the domain
air_elements = np.arange(0, 100)
xfm_elements = np.arange(100, 120)
foam_elements = np.arange(120, num_nodes)
subdomains = {air: air_elements, xfm: xfm_elements, foam:foam_elements}
check_material_compability(subdomains)


# initialize the assembler
assembler = Assembler(dof_handler, bases, subdomains, dtype=np.complex128)

K_g= assembler.assemble_material_K(omega)  # global stiffness matrix with material attribution
M_g= assembler.assemble_material_M(omega)  # global mass matrix with material attribution
# print("K_g:", assembler.get_matrix_in_array(K_g))

# construct linear system
left_hand_matrix = K_g-M_g
# plot_matrix_partten(left_hand_matrix)

# print(display_matrix_in_array(left_hand_matrix))
#  natural boundary condition   
nature_bcs = {'type': 'velocity', 'value': 1*np.exp(-1j*omega), 'position': 0}
right_hand_vector = assembler.assemble_nature_bc(nature_bcs)
# print(right_hand_vector)

# solver the linear system
linear_solver = LinearSolver(dof_handler)
linear_solver.solve(left_hand_matrix, right_hand_vector)
sol = linear_solver.u
# print("solution:", abs(sol))


# plot the solution
# post_process = PostProcessField(mesh.nodes, r'1D Helmholtz (2000$Hz$)')
# post_process.plot_sol((np.real(sol), f'FEM ($p=3$)', 'solid'))
# post_process.display_layers(nodes[99], nodes[119])
# plt.show()

# plot the solution as animation as with change of time
fig = plt.figure()
ax = plt.axes(xlim=(-1, 1), ylim=(-300, 300))
line, = ax.plot([], [], lw=2)

time = np.linspace(0, 1, 201)

plt.ion()
for i in range(len(time)):
    plt.clf()
    plt.plot(nodes, np.real(sol*np.exp(-1j*omega*i)))
    plt.pause(0.1)

