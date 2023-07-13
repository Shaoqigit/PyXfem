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
    def read_mesh(self):
        pass


class Mesh1D(BaseMesh):

    def __init__(self, nodes, elem_connect):
        self.nodes = nodes
        self.elem_connect = elem_connect
        self.mesh = self.read_mesh()

    def read_mesh(self):
        """read mesh file"""
        elems = {}
        for i in range(len(self.elem_connect)):
            elems[i] = np.array([self.nodes[i], self.nodes[i+1]])
        return elems

    def get_num_nodes(self):
        """return number of nodes"""
        return len(self.nodes)
    
    def get_num_elems(self):
        """return number of elements"""
        return len(self.elem_connect)
    
    @property
    def connectivity(self):
        """return connectivity"""
        return self.elem_connect
    


    


            