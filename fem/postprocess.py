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

# simple postprocessor for plotting and error computation

import numpy as np
import matplotlib.pyplot as plt


class PostProcessField(object):

    def __init__(self, x_nodes, title):
        self.x_nodes = x_nodes
        self.title = title
        self.set_figure('Position(m)', 'Pressure(Pa)')

    def set_figure(self, xaxis, yaxis):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.set_title(self.title)
        self.ax.set_xlabel(xaxis)
        self.ax.set_ylabel(yaxis)
        
    
    def plot_sol(self, *sols):
        for sol in sols:
            self.ax.plot(self.x_nodes, sol[0], label=sol[1], linestyle=sol[2])
        
        self.ax.legend()

    def display_layers(self, *layers_pos):
        for pos in layers_pos:
            self.ax.axvline(x=pos, ls='--', c='k')


    def compute_error(self, sol, ana_sol, remove=None):
        # relative error

        # Compute the differences between the solutions
        differences = np.abs(ana_sol - sol)

        # Square the differences
        squared_differences = differences ** 2

        # Sum up the squared differences
        sum_squared_differences = np.sum(squared_differences)

        # Divide by the number of nodes
        num_nodes = ana_sol.size
        mean_squared_difference = sum_squared_differences / num_nodes

        # Compute analytical
        abs_analytical = np.sum(np.abs(ana_sol)**2)/num_nodes

        # Relative error
        l2_error = np.sqrt(mean_squared_difference/abs_analytical)

        return l2_error
    
    # def compute_L2_error(self, mesh, sol, ana_sol):
    #     # L2 error
    #     # rough error computation
    #     error = 0
    #     for i in range(len(mesh.elements)):
    #         error += 
    #     return error