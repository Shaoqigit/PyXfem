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


class Mesh1D(BaseMesh):

    def __init__(self, nodes, elem_connect):
        self.nodes = nodes
        self.elem_connect = elem_connect

    def get_mesh(self):
        """dict of element number and nodes coordinates"""
        elems = {}
        for i in range(len(self.elem_connect)):
            elems[i] = np.array([self.nodes[i], self.nodes[i+1]])
        return elems
    
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
        

    def get_nb_nodes(self):
        """return number of nodes"""
        return len(self.nodes)
    
    def get_nb_elems(self):
        """return number of elements"""
        return len(self.elem_connect)
    
    def get_nodes_from_elem(self, elem):
        """return nodes of corresoinding element"""
        return self.get_mesh()[elem]
    
    


    


            