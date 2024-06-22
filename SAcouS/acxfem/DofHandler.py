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

# dofhandle.py
# deal with mesh and basis, return global dof index t

import numpy as np
from collections import defaultdict

from .Mesh import Mesh1D


class DofHandler1D:

  def __init__(self, mesh, bases) -> None:
    if isinstance(mesh, Mesh1D):
      self.mesh = mesh
      self.bases = bases
      self.connect = mesh.connectivity
    else:
      raise TypeError("mesh must be 1D mesh")

  def get_nb_dofs(self):
    """return global dof"""
    nb_dofs = 0
    for basis in self.bases:
      if basis.is_discontinue:
        num_discontiue = self.basis.interface
        return self.mesh.get_nb_elems() * self.basis.get_order(
        ) + 1 + num_discontiue * (self.basis.get_order() + 1)
      else:
        nb_dofs += basis.get_order()
    return nb_dofs + 1

  def get_global_dofs(self):
    """return local dof
        return: global dof index
        [node_1_index, node_2_index, internal_dof_index],
        [.....],
        [.....]]"""
    global_dof = []
    internal_dof_index_start = self.mesh.get_nb_nodes()
    for i, basis in enumerate(self.bases):
      # print("internal_dof_index_start", internal_dof_index_start)
      if basis.nb_internal_dofs == 0:
        global_dof.append(self.mesh.connectivity[i])
        continue
      elem = self.mesh.connectivity[i]
      global_dof.append(
          np.hstack((np.array([elem[0], elem[1]]),
                     np.array([
                         internal_dof_index_start + j
                         for j in range(0,
                                        basis.get_order() - 1)
                     ]))))
      internal_dof_index_start += basis.get_order() - 1

    return global_dof

  @property
  def nb_external_dofs(self):
    return self.mesh.get_nb_nodes()

  @property
  def nb_internal_dofs(self):
    return self.get_nb_dofs() - self.nb_external_dofs


class DofHandler1DMutipleVariable(DofHandler1D):

  def __init__(self, mesh, *mbases) -> None:
    if isinstance(mesh, Mesh1D):
      self.mesh = mesh
      self.nb_var = len(mbases)
      self.whole_bases = []
      self.var_name = []
      for bases in mbases:
        self.whole_bases += bases
        if bases[0].label not in self.var_name:
          self.var_name.append(bases[0].label)
      self.connect = mesh.connectivity

    else:
      raise TypeError("mesh must be 1D mesh")

  def get_nb_dofs(self):
    """return global dof"""
    nb_dofs = 0
    for basis in self.whole_bases:
      if basis.is_discontinue:
        num_discontiue = self.basis.interface
        return self.mesh.get_nb_elems() * self.basis.get_order(
        ) + 1 + num_discontiue * (self.basis.get_order() + 1)
      else:
        nb_dofs += basis.get_order()
    return (nb_dofs + self.nb_var)

  def get_global_dofs(self):
    """return local dof
        return: global dof index
        [node_1_index, node_2_index, internal_dof_index],
        [.....],
        [.....]]"""
    global_dof = []
    internal_dof_index_start = self.mesh.get_nb_nodes() * self.nb_var

    new_connect = self.connect
    whole_connect = self.connect
    for i in range(1, self.nb_var):
      index_1 = new_connect[:, 1] + self.connect.shape[0] * i
      index_2 = new_connect[:, 1] + 1 + self.connect.shape[0] * i
      new_connect = np.vstack((index_1, index_2)).T
      whole_connect = np.vstack((whole_connect, new_connect))

    for i, basis in enumerate(self.whole_bases):
      if basis.nb_internal_dofs == 0:
        global_dof.append(whole_connect[i])
        continue
      elem = whole_connect[i]
      global_dof.append(
          np.hstack((np.array([elem[0], elem[1]]),
                     np.array([
                         internal_dof_index_start + j
                         for j in range(0,
                                        basis.get_order() - 1)
                     ]))))
      internal_dof_index_start += basis.get_order() - 1

    return global_dof

  def base4global_dofs(self):
    global_dofs = self.get_global_dofs()
    base4dofs = {}
    for i, bases in enumerate(self.whole_bases):
      index = base4dofs.get(bases.label, [])
      index.append(global_dofs[i])
      base4dofs[bases.label] = index
    return base4dofs

  def get_global_dofs_by_base(self, label):
    base4dofs = self.base4global_dofs()
    return base4dofs.get(label, None)

  def mesh2dof(self, position, var):
    index_var = self.var_name.index(var)
    for i, node in enumerate(self.mesh.nodes):
      if node == position:
        return i + index_var * self.mesh.get_nb_nodes()

  @property
  def nb_external_dofs(self):
    return self.mesh.get_nb_nodes() * self.nb_var

  @property
  def nb_internal_dofs(self):
    return self.get_nb_dofs() - self.nb_external_dofs


class GeneralDofHandler1D(DofHandler1D):

  def __init__(self, var_name, *mbases) -> None:
    self.nb_var = len(mbases)
    self.var_name = var_name
    self.whole_bases = []
    # coupling bases must have same dimension
    if len(mbases) > 1:
      for bases in mbases:
        if len(bases) != len(mbases[0]):
          raise ValueError("the dimensionn of coupling bases must be the same")

    for bases in mbases:
      self.whole_bases += bases
    self.nb_elem = len(mbases[0])
    self.num_nodes = self.nb_elem + 1

  def get_nb_dofs(self):
    """return global dof"""
    nb_dofs = 0
    for basis in self.whole_bases:
      if basis.is_discontinue:
        num_discontiue = self.basis.interface
        return self.num_nodes * self.basis.get_order() + 1 + num_discontiue * (
            self.basis.get_order() + 1)
      else:
        nb_dofs += basis.get_order()
    return (nb_dofs + self.nb_var)

  def get_global_dofs(self):
    """return local dof
        return: global dof index
        [node_1_index, node_2_index, internal_dof_index],
        [.....],
        [.....]]"""
    global_dof = []
    internal_dof_index_start = self.num_nodes * self.nb_var
    elem_connec1 = np.arange(0, self.nb_elem)
    elem_connec2 = np.arange(1, self.num_nodes)
    connect = np.vstack((elem_connec1, elem_connec2)).T
    new_connect = connect
    whole_connect = connect
    for i in range(1, self.nb_var):
      index_1 = new_connect[:, 1] + connect.shape[0] * i
      index_2 = new_connect[:, 1] + 1 + connect.shape[0] * i
      new_connect = np.vstack((index_1, index_2)).T
      whole_connect = np.vstack((whole_connect, new_connect))

    for i, basis in enumerate(self.whole_bases):
      if basis.nb_internal_dofs == 0:
        global_dof.append(whole_connect[i])
        continue
      elem = whole_connect[i]
      global_dof.append(
          np.hstack((np.array([elem[0], elem[1]]),
                     np.array([
                         internal_dof_index_start + j
                         for j in range(0,
                                        basis.get_order() - 1)
                     ]))))
      internal_dof_index_start += basis.get_order() - 1

    return global_dof

  def base4global_dofs(self):
    global_dofs = self.get_global_dofs()
    base4dofs = {}
    for i, bases in enumerate(self.whole_bases):
      index = base4dofs.get(bases.label, [])
      index.append(global_dofs[i])
      base4dofs[bases.label] = index
    return base4dofs

  def get_global_dofs_by_base(self, label):
    base4dofs = self.base4global_dofs()
    return base4dofs.get(label, None)

  def mesh2dof(self, position, var):
    index_var = self.var_name.index(var)
    for i, node in enumerate(self.mesh.nodes):
      if node == position:
        return i + index_var * self.mesh.get_nb_nodes()

  @property
  def num_external_dofs(self):
    return self.num_nodes * self.nb_var

  @property
  def num_internal_dofs(self):
    return self.get_nb_dofs() - self.num_external_dofs


class FESpace:

  def __init__(self, mesh, subdomains, *all_bases) -> None:
    self.mesh = mesh
    self.subdomains = subdomains
    self.nb_var = len(all_bases)

    self.whole_bases = [base for bases in all_bases for base in bases]
    self.var_names = [base.label for base in self.whole_bases]

    self.elem_mat = {}
    self.subdoamin_start_index = {}
    for mat, elems in self.subdomains.items():
      self.subdoamin_start_index[mat] = elems[0]
      self.elem_mat.update({elem: mat for elem in elems})

  @property
  def nb_external_dofs(self):
    return self.mesh.get_nb_nodes() * self.nb_var

  @property
  def nb_internal_dofs(self):
    return self.get_nb_dofs() - self.nb_external_dofs

  def get_nb_dofs(self):
    # to be modified to adapet to 2D/3D
    nb_nodes = self.mesh.get_nb_nodes()
    nb_dofs = nb_nodes * self.nb_var
    for basis in self.whole_bases:
      if basis.is_discontinue:
        num_discontiue = self.basis.interface
        return nb_nodes * self.basis.get_order() + 1 + num_discontiue * (
            self.basis.get_order() + 1)
      else:
        nb_dofs += basis.nb_internal_dofs
    return nb_dofs

  def get_global_dofs(self):
    """return local dof
        return: global dof index
        [node_1_index, node_2_index, internal_dof_index],
        [.....],
        [.....]]"""
    global_dof = []
    internal_dof_index_start = self.mesh.get_nb_nodes() * self.nb_var
    connect = self.mesh.connectivity
    new_connect = connect
    whole_connect = connect
    for i in range(1, self.nb_var):
      index_1 = new_connect[:, 1] + connect.shape[0] * i
      index_2 = new_connect[:, 1] + 1 + connect.shape[0] * i
      new_connect = np.vstack((index_1, index_2)).T
      whole_connect = np.vstack((whole_connect, new_connect))

    for i, basis in enumerate(self.whole_bases):
      if basis.nb_internal_dofs == 0:
        global_dof.append(whole_connect[i])
        continue
      elem = whole_connect[i]
      global_dof.append(
          np.hstack((np.array([elem[0], elem[1]]),
                     np.array([
                         internal_dof_index_start + j
                         for j in range(0,
                                        basis.get_order() - 1)
                     ]))))
      internal_dof_index_start += basis.get_order() - 1

    return global_dof

  def base4global_dofs(self):
    """
    return the global dofs for each base
    """
    global_dofs = self.get_global_dofs()
    base2dofs = defaultdict(list)
    for i, bases in enumerate(self.whole_bases):
      base2dofs[bases.label].append(global_dofs[i])
    return base2dofs

  def get_global_dofs_by_base(self, label):
    base4dofs = self.base4global_dofs()
    return base4dofs.get(label, None)

  def get_dofs_from_var_coord(self, coord, var):
    self.mat2dofs = defaultdict(
        dict)    # {mat: {coord: dofs}} for BC imposition
    for mat, elems in self.subdomains.items():
      start_index = self.subdoamin_start_index[mat]
      for elem in elems:
        local_dofs = self.mesh.connectivity[elem]
        local_coord = [self.mesh.num_node2coord[node] for node in local_dofs]
        if mat.TYPE in ['Air', 'Fluid']:
          self.mat2dofs['Pf'].update(dict(zip(local_coord,
                                              local_dofs)))    # neat opeartor
        elif mat.TYPE in ['Poroelastic']:
          self.mat2dofs['Pb'].update(dict(zip(local_coord, local_dofs)))
          local_dofs_2 = local_dofs + self.mesh.get_nb_nodes() - start_index
          self.mat2dofs['Ux'].update(dict(zip(local_coord, local_dofs_2)))
    return self.mat2dofs[var][coord]

  def mesh2dof(self, position, var):
    index_var = self.var_name.index(var)
    for i, node in enumerate(self.mesh.nodes):
      if node == position:
        return i + index_var * self.mesh.get_nb_nodes()
