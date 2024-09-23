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
from typing import Union

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
  def get_mesh_coordinates(self):
    """return element to nodes coordinates"""
    pass

  @abstractmethod
  def refine_mesh(self, times):
    pass

  @abstractmethod
  def plotmesh(self):
    pass

  @property
  def connectivity(self):
    """return connectivity"""
    return self.elem_connect

  def node2elem(self, node):
    """return element number from node number"""
    for i in range(len(self.elem_connect)):
      if node in self.elem_connect[i]:
        return i
    raise ValueError("node not in mesh")

  def get_nb_nodes(self):
    """return number of nodes"""
    return len(self.nodes)

  def get_nb_elems(self):
    """return number of elements"""
    return len(self.elem_connect)

  def get_nodes_from_elem(self, elem):
    """return nodes of corresoinding element"""
    return self.get_mesh_coordinates()[elem]

  @property
  def num_node2coord(self):
    num_node2coord2 = {}
    """return node number from coordinate"""
    for i, coord in enumerate(self.nodes):
      if isinstance(coord, np.ndarray):    # for 2/3D mesh
        coord = tuple(coord)
      num_node2coord2[i] = coord
    return num_node2coord2


class Mesh1D(BaseMesh):

  def __init__(self, nodes, elem_connect):
    self.nodes = nodes
    self.nb_nodes = len(nodes)
    self.elem_connect = elem_connect
    self.dim = 1
    self.node_index = np.arange(self.nb_nodes)
    self._subdomains = None

  def set_subdomains(self, subdomains):
    self._subdomains = subdomains

  @property
  def subdomains(self):
    return self._subdomains

  def get_mesh_coordinates(self):
    """dict of element number and nodes coordinates"""
    elems = {}
    for i in range(len(self.elem_connect)):
      elems[i] = np.array([self.nodes[i], self.nodes[i + 1]])
    return elems

  def get_min_size(self):
    """return minimum size of elements"""
    return min(np.diff(self.nodes))

  @property
  def coord2node_num(self):
    """return element number from node number"""
    coord2node_num = {}
    for i, coord in enumerate(self.nodes):
      coord2node_num[coord] = i
    return coord2node_num

  def refine_mesh(self, times):
    """refine mesh"""
    for _ in range(times):
      new_nodes = []
      new_elem_connect = []
      for i in range(len(self.elem_connect)):
        new_nodes.append(self.nodes[i])
        new_nodes.append(0.5 * (self.nodes[i] + self.nodes[i + 1]))
        new_elem_connect.append(np.array([2 * i, 2 * i + 1]))
        new_elem_connect.append(np.array([2 * i + 1, 2 * i + 2]))
      new_nodes.append(self.nodes[-1])
      self.nodes = np.array(new_nodes)
      self.elem_connect = np.array(new_elem_connect)
      self.nb_nodes = len(self.nodes)

  def plotmesh(self, withnode=False, withnodeid=False):
    """
    plot the 1d mesh

    Parameters
    ----------
    withnode : boolean
        True to show the node
    withnodeid : boolean
        True to show the node id
    """
    y = np.zeros(len(self.nodes))
    plt.figure()
    plt.plot(self.nodes, y, 'k')
    if withnode:
      plt.plot(self.nodes, y, 'r*')
    if withnodeid:
      for i, node in enumerate(self.nodes):
        x = self.nodes[i]
        plt.text(x, 0.0, '%d' % (i + 1))
    plt.xlabel('X', fontsize=14)
    plt.ylabel('Y', fontsize=14)
    plt.show()


import meshio


class Mesh2D(BaseMesh):

  def __init__(self, nodes, elem_connect, edge_connect=None):
    self.nodes = nodes
    self.elem_connect = elem_connect
    self.nb_nodes = len(self.nodes)
    self.node_index = np.arange(self.nb_nodes)
    self.exterior_facets = edge_connect    # [node1, node2]
    self.mesh_center = np.mean(self.nodes, axis=0)    # center of the mesh
    self.dim = 2
    self._subdomains = None

  def set_subdomains(self, subdomains):
    self._subdomains = subdomains

  @property
  def subdomains(self):
    return self._subdomains

  def plotmesh(self, withnode=False, withnodeid=False, withedgeid=False):
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
    if withedgeid:
      for i, edge in enumerate(self.exterior_facets):
        x, y = self.nodes[edge[0]]
        x1, y1 = self.nodes[edge[1]]
        plt.text(0.5 * (x + x1), 0.5 * (y + y1), '%d' % (i + 1))
    plt.xlabel('X', fontsize=14)
    plt.ylabel('Y', fontsize=14)
    plt.show()

  def get_mesh_coordinates(self):
    """dict of element number and nodes coordinates"""
    return {i: self.nodes[conn] for i, conn in enumerate(self.elem_connect)}

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
        new_nodes.extend([
            self.nodes[n1], n12, n31, n12, self.nodes[n2], n23, n31, n23,
            self.nodes[n3]
        ])
        new_elem_connect.append(
            [len(new_nodes) - 9,
             len(new_nodes) - 8,
             len(new_nodes) - 7])
        new_elem_connect.append(
            [len(new_nodes) - 6,
             len(new_nodes) - 5,
             len(new_nodes) - 4])
        new_elem_connect.append(
            [len(new_nodes) - 3,
             len(new_nodes) - 2,
             len(new_nodes) - 1])
        self.nodes = np.array(new_nodes)
        self.elem_connect = np.array(new_elem_connect)
        self.nb_nodes = len(self.nodes)

  def compute_normal(self, edge):
    """compute normal of the element"""
    assert (len(edge) == 2)
    edge = edge
    node_1_coord = self.nodes[edge[0]]    # in physical space
    node_2_coord = self.nodes[edge[1]]
    tangent = np.array(
        [node_2_coord[0] - node_1_coord[0], node_2_coord[1] - node_1_coord[1]])
    normal = np.array([tangent[1], -tangent[0]])
    normal = normal / np.linalg.norm(normal)
    return normal


class Mesh3D(Mesh2D):

  def __init__(self, nodes, elem_connect, surface_connect):
    super().__init__(nodes, elem_connect, surface_connect)
    self.dim = 3
    self.surface_connect = surface_connect
    self.mesh_center = np.mean(self.nodes, axis=0)

  def get_mesh_coordinates(self):
    return super().get_mesh_coordinates()

  def plotmesh(self, withnode=False, withnodeid=False, withedgeid=False):
    raise NotImplementedError("3D mesh plot not implemented yet")

  def compute_normal(self, facet):
    """compute normal of the element"""
    assert (len(facet) == 3)
    facet = facet
    node_1_coord = self.nodes[facet[0]]    # in physical space
    node_2_coord = self.nodes[facet[1]]
    node_3_coord = self.nodes[facet[2]]
    vect_1 = node_2_coord - node_1_coord
    vect_2 = node_3_coord - node_1_coord
    # compute the normal vector with cross product
    normal = np.cross(vect_1, vect_2)
    normal = normal / np.linalg.norm(normal)
    element_center = np.mean([node_1_coord, node_2_coord, node_3_coord],
                             axis=0)
    outwards = element_center - self.mesh_center
    if np.dot(normal, outwards) < 0:
      normal = -normal
    return normal


def mesh_constructor(dim, nodes, elem_connect, edge_connect=None):
  if dim == 1:
    return Mesh1D(nodes, elem_connect)
  elif dim == 2:
    return Mesh2D(nodes, elem_connect, edge_connect)
  elif dim == 3:
    return Mesh3D(nodes, elem_connect, edge_connect)
  else:
    raise ValueError("dimension not supported")


GMSH_DIM_FACET_MAP = {2: 'triangle', 3: 'tetra'}
GMSH_DIM_EDGE_MAP = {2: 'line', 3: 'triangle'}


class MeshReader:

  def __init__(self, mesh_file_name, dim=2):
    self.extension = mesh_file_name.split('.')[-1]
    self.dim = dim
    self.meshio_object = meshio.read(mesh_file_name)

  def get_mesh(self) -> BaseMesh:
    if self.extension == 'msh':
      # version 2.2 without saving all parameters
      nodes = self.meshio_object.points[:, :self.dim]
      for elem in self.meshio_object.cells:
        if elem.type == GMSH_DIM_FACET_MAP[self.dim]:
          elem_connect = elem.data
        elif elem.type == GMSH_DIM_EDGE_MAP[self.dim]:
          facet_connect = elem.data
      return mesh_constructor(self.dim, nodes, elem_connect, facet_connect)

  def init_subdomains(self, mesh: BaseMesh,
                      physical_tag2material: dict) -> dict:
    """return subdomains"""
    subdomains = {}
    for tag, material in physical_tag2material.items():
      subdomains[material] = self.get_elem_by_physical(tag)
    mesh.set_subdomains(subdomains)
    return subdomains

  def get_elem_by_physical(self, physical_tag: Union[str, int]) -> np.ndarray:
    """return element number by physical tag
    3D elements in dim=3
    2D elements in dim=2
    """
    if isinstance(physical_tag, str):
      elem_tag = int(self.meshio_object.field_data[physical_tag][0])
    else:
      elem_tag = physical_tag
    elem_index = np.where(self.meshio_object.cell_data_dict['gmsh:physical'][
        GMSH_DIM_FACET_MAP[self.dim]] == elem_tag)
    return elem_index[0]

  def get_facet_by_physical(self, physical_tag: Union[str, int]) -> np.ndarray:
    """return facet number by physical tag
    2D facets in dim=3
    1D facets in dim=2
    """
    if isinstance(physical_tag, str):
      edge_tag = int(self.meshio_object.field_data[physical_tag][0])
    else:
      edge_tag = physical_tag
    edge_index = np.where(self.meshio_object.cell_data_dict['gmsh:physical'][
        GMSH_DIM_EDGE_MAP[self.dim]] == edge_tag)
    return edge_index[0]

  def get_vertices_by_physical(self, physical_tag: Union[str,
                                                         int]) -> np.ndarray:
    """return vertices number by physical tag"""
    raise NotImplementedError("vertices not implemented yet")


if __name__ == "__main__":
  mesh = Mesh2D()
  mesh_square = mesh.read_mesh("../../tests/mesh/square_1.msh")
  mesh.plotmesh(withnode=True, withnodeid=True)
