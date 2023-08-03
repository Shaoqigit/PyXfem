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
        # rough error computation
        sol = sol[:remove]
        ana_sol = ana_sol[:remove]
        error = np.mean(np.abs(sol-ana_sol)/np.abs(ana_sol))
        # error = np.mean(np.abs(sol-ana_sol))
        return error