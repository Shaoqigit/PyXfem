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

# mesh.py: generate mesh data dictionsary and refine mesh function


import numpy as np
import matplotlib.pyplot as plt
from abc import ABCMeta, abstractmethod

class BaseMesh(metaclass=ABCMeta):
    """base abstract mesh class"""

    def parser_mesh(mesh_file):
        """parse mesh file"""
        with open(mesh_file, 'r') as f:
            lines = f.readlines()
        return lines
   
    @abstractmethod
    def get_mesh(self):
        pass

    @abstractmethod
    def refine_mesh(self, times):
        pass

    @abstractmethod
    def plotmesh(self):
        pass


class Mesh1D(BaseMesh):

    def __init__(self, nodes, elem_connect):
        self.nodes = nodes
        self.nb_nodes = len(nodes)
        self.elem_connect = elem_connect

    def get_mesh(self):
        """dict of element number and nodes coordinates"""
        elems = {}
        for i in range(len(self.elem_connect)):
            elems[i] = np.array([self.nodes[i], self.nodes[i+1]])
        return elems
    
    def get_min_size(self):
        """return minimum size of elements"""
        return min(np.diff(self.nodes))
    
    def node2elem(self, node):
        """return element number from node number"""
        for i in range(len(self.elem_connect)):
            if node in self.elem_connect[i]:
                return i
        raise ValueError("node not in mesh")

    @property
    def num_node2coord(self):
        num_node2coord2 = {}
        """return node number from coordinate"""
        for i, coord in enumerate(self.nodes):
            num_node2coord2[i] = coord
        return num_node2coord2
    
    @property
    def coord2node_num(self):
        """return element number from node number"""
        coord2node_num = {}
        for i, coord in enumerate(self.nodes):
            coord2node_num[coord] = i
        return coord2node_num
    

    @property
    def connectivity(self):
        """return connectivity"""
        return self.elem_connect
    

    def refine_mesh(self, times):
        """refine mesh"""
        for _ in range(times):
            new_nodes = []
            new_elem_connect = []
            for i in range(len(self.elem_connect)):
                new_nodes.append(self.nodes[i])
                new_nodes.append(0.5*(self.nodes[i] + self.nodes[i+1]))
                new_elem_connect.append(np.array([2*i, 2*i+1]))
                new_elem_connect.append(np.array([2*i+1, 2*i+2]))
            new_nodes.append(self.nodes[-1])
            self.nodes = np.array(new_nodes)
            self.elem_connect = np.array(new_elem_connect)
            self.nb_nodes = len(self.nodes)

    def plotmesh(self,withnode=False,withnodeid=False):
        """
        plot the 1d mesh

        Parameters
        ----------
        withnode : boolean
            True to show the node
        withnodeid : boolean
            True to show the node id
        """
        y=np.zeros(len(self.nodes))
        plt.figure()
        plt.plot(self.nodes,y,'k')
        if withnode:
            plt.plot(self.nodes,y,'r*')
        if withnodeid:
            for i, node in enumerate(self.nodes):
                x=self.nodes[i]
                plt.text(x,0.0,'%d'%(i+1))
        plt.xlabel('X',fontsize=14)
        plt.ylabel('Y',fontsize=14)
        plt.show()
        

    def get_nb_nodes(self):
        """return number of nodes"""
        return len(self.nodes)
    
    def get_nb_elems(self):
        """return number of elements"""
        return len(self.elem_connect)
    
    def get_nodes_from_elem(self, elem):
        """return nodes of corresoinding element"""
        return self.get_mesh()[elem]

import meshio
class Mesh2D(BaseMesh):
    def __init__(self):
        self.nodes = []
        self.elem_connect = []
        self.nb_nodes = 0

    def read_mesh(self, mesh_file):
        mesh = meshio.read(mesh_file)
        self.nodes = mesh.points[:,:2]
        self.elem_connect = mesh.cells[2].data
        self.nb_elmes = len(self.elem_connect)
        self.nb_nodes = len(self.nodes)

    def plotmesh(self, withnode=False, withnodeid=False):
        """
        plot the 2d mesh

        Parameters
        ----------
        withnode : boolean
            True to show the node
        withnodeid : boolean
            True to show the node id
        """
        plt.figure()
        for elem in self.elem_connect:
            x = self.nodes[elem][:, 0]
            y = self.nodes[elem][:, 1]
            plt.plot(np.append(x, x[0]), np.append(y, y[0]), 'k')
        if withnode:
            plt.plot(self.nodes[:, 0], self.nodes[:, 1], 'r*')
        if withnodeid:
            for i, node in enumerate(self.nodes):
                x, y = self.nodes[i]
                plt.text(x, y, '%d' % (i + 1))
        plt.xlabel('X', fontsize=14)
        plt.ylabel('Y', fontsize=14)
        plt.show()


    def get_mesh(self):
        """dict of element number and nodes coordinates"""
        elems = {}
        for i in range(len(self.elem_connect)):
            elems[i] = self.nodes[self.elem_connect[i]]
        return elems

    def refine_mesh(self, times):
        # refine the 2D mesh
        for _ in range(times):
            new_nodes = []
            new_elem_connect = []
            for elem in self.elem_connect:
                n1, n2, n3 = elem
                n12 = 0.5 * (self.nodes[n1] + self.nodes[n2])
                n23 = 0.5 * (self.nodes[n2] + self.nodes[n3])
                n31 = 0.5 * (self.nodes[n3] + self.nodes[n1])
                new_nodes.extend([self.nodes[n1], n12, n31, n12, self.nodes[n2], n23, n31, n23, self.nodes[n3]])
                new_elem_connect.append([len(new_nodes)-9, len(new_nodes)-8, len(new_nodes)-7])
                new_elem_connect.append([len(new_nodes)-6, len(new_nodes)-5, len(new_nodes)-4])
                new_elem_connect.append([len(new_nodes)-3, len(new_nodes)-2, len(new_nodes)-1])
                self.nodes = np.array(new_nodes)
                self.elem_connect = np.array(new_elem_connect)
                self.nb_nodes = len(self.nodes)
        

    
if __name__ == "__main__":
    mesh = Mesh2D()
    mesh_square = mesh.read_mesh("../../tests/mesh/square_1.msh")
    mesh.plotmesh(withnode=True, withnodeid=True)
    
    


    


            