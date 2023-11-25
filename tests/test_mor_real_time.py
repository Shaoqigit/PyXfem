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

from acxfem.basis import Lobbato1DElement
from acxfem.mesh import Mesh1D
from acxfem.dofhandler import DofHandler1D
from acxfem.assembly import Assembler
from acxfem.materials import Air, Fluid, EquivalentFluid
from acxfem.utilities import check_material_compability, display_matrix_in_array, plot_matrix_partten
from acxfem.solver import LinearSolver

from acxmor.modal_reduction import EigenSolver, ModalReduction
from acxfem.postprocess import PostProcessField, PostProcessFRF
from analytical.fluid_sol import ImpedenceKundltTube


def FEM_model(omega, assembler, nature_bcs):
    assembler.initial_matrix()
    K_g= assembler.assemble_material_K(omega)  # global stiffness matrix with material attribution
    M_g= assembler.assemble_material_M(omega)  # global mass matrix with material attribution
    right_hand_vector = assembler.assemble_nature_bc(nature_bcs)
    left_hand_matrix = K_g-M_g

    # solver the linear system
    linear_solver = LinearSolver(dof_handler)
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
    # Harmonic Acoustic problem define the frequency
    # Set the mean frequency and standard deviation
    mean_frequency = 1000  # in Hz
    std_deviation = 10

    # Generate random frequencies using a normal distribution
    num_frequencies = 500
    freq = np.random.normal(loc=mean_frequency, scale=std_deviation, size=num_frequencies)
    omegas = 2 * np.pi * freq  # angular frequency

    num_elem = 500  # number of elements
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
    nature_bcs = {'type': 'total_displacement', 'value': 1, 'position': 0}
    impedence_bcs = {'type': 'impedence', 'value': 0., 'position': num_elem}

    import time
    # start = time.time()
    # sol_fem = np.zeros((100,num_nodes), dtype=np.complex128)
    # for i, omega in enumerate(omegas):
    #     sol_fem[i,:] = FEM_model(omega, assembler, nature_bcs)[int(num_nodes/2)]
    # end = time.time()
    # print("FEM solving time: ", end-start)

    # ====================== Modal Reduction ======================
    start = time.time()
    nb_modes = 50
    assembler.initial_matrix()
    # import pdb; pdb.set_trace()
    K_w= assembler.assemble_material_K(omega=1)  # global stiffness matrix no frequency dependent material
    M_w= assembler.assemble_material_M(omega=1)  # global mass matrix: the eigen freq should be omega^2
    # construct modal domain
    eigen_value_solver = EigenSolver(dof_handler)
    eig_omega_sq, modes = eigen_value_solver.solve(K_w, M_w, nb_modes)
    # eig_freqs = np.sqrt(eig_omega_sq)/(2*np.pi)
    
    modal_reduction_method = ModalReduction(K_w, M_w, modes)

    #  natural boundary condition   
    right_hand_vector = assembler.assemble_nature_bc(nature_bcs)
    # ====================== Reduced System ======================
    # reduced the system
    K_r = modal_reduction_method.projection(K_w)
    M_r = modal_reduction_method.projection(M_w)
    f_r = modal_reduction_method.projection(right_hand_vector)

    sol_modal = np.zeros((num_frequencies,num_nodes), dtype=np.complex128)
    # for i, omega in enumerate(omegas):
    #     sol_modal[i, :] = reduction_model(omega, K_r, M_r, f_r)
    # end = time.time()
    # print("Modal reduction solving time: ", end-start)
        # ====================== Analytical Solution ======================

    # plot the solution
    x = nodes
    y = sol_modal[0]

    # You probably won't need this if you're embedding things in a tkinter plot...
    plt.ion()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1, = ax.plot(x, y, 'r-') # Returns a tuple of line objects, thus the comma
    ax.set_title('Real-Time Kundlt Tube Pressure (FEM)')
    ax.set_xlabel('Position(m)')
    ax.set_ylim(-0.5, 0.5)
    ax.set_ylabel('Pressure(Pa)')
    for i, omega in enumerate(omegas):
        print ("frequency: ", omega/(2*np.pi))
        line1.set_ydata(reduction_model(omega, K_r, M_r, f_r))
        fig.canvas.draw()
        fig.canvas.flush_events()
        time.sleep(0.1)


    # import sounddevice as sd

    # def plot_spectrum(data, sample_rate):
    #     # Perform FFT
    #     fft_result = np.fft.fft(data)
    #     # Calculate the corresponding frequencies
    #     frequencies = np.fft.fftfreq(len(data), d=1/sample_rate)
    #     # Plot the magnitude spectrum
    #     plt.plot(frequencies, np.abs(fft_result))
    #     plt.title('Frequency Spectrum')
    #     plt.xlabel('Frequency (Hz)')
    #     plt.ylabel('Magnitude')
    #     plt.show()

    # def callback(indata, frames, time, status):
    #     if status:
    #         print(status)
    #     # Plot the frequency spectrum for each callback
    #     plot_spectrum(indata[:, 0], sample_rate)

    # # Set the sample rate and duration
    # sample_rate = 44100  # You can adjust this based on your needs
    # duration = 5  # Set the duration of the recording in seconds

    # # Capture microphone input
    # with sd.InputStream(callback=callback, channels=1, samplerate=sample_rate):
    #     sd.sleep(int(duration * 1000))
    # plt.show(block=True)