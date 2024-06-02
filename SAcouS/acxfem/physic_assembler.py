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

# assembly the global/partial matrices according to the physic of the components
from SAcouS.acxfem.mesh import Mesh1D
from SAcouS.acxfem.dofhandler import DofHandler1D

import numpy as np
from numba import jit
from scipy.sparse import csr_array, coo_matrix, lil_array
from multiprocessing import Pool


def get_indeces(*dofs):
  if len(dofs) == 1:
    return np.transpose(np.meshgrid(dofs[0], dofs[0])).reshape(-1, 2)
  elif len(dofs) == 2:
    return np.transpose(np.meshgrid(dofs[0], dofs[1])).reshape(-1, 2)
  else:
    raise ValueError("wrong number of arguments")


# def assembly_global_matrix(bases, elem_mat, omega, var=None)
class BaseAssembler:

  def __init__(self, fe_space, omega, dtype) -> None:
    """
        General assembler for Helmholtz equation
        bases: list of basis
        subdomains: dict of subdomains
        dtype: data type of linear system"""
    self.fe_space = fe_space
    self.nb_global_dofs = fe_space.get_nb_dofs()
    self.dtype = dtype
    self.elem_mat = fe_space.elem_mat
    self.omega = omega

  def initial_matrix(self):
    self.K = 0.
    self.M = 0.

  def assemble_material_K(self, bases, var=None):
    if var is None:
      dofs_index = self.fe_space.get_global_dofs()
    else:
      dofs_index = self.fe_space.get_global_dofs_by_base(var)

    rows = []
    cols = []
    data = []
    for i, (dofs, basis) in enumerate(zip(dofs_index, bases)):
      local_indices = np.stack(
          (np.repeat(basis.local_dofs_index, len(basis.local_dofs_index)),
           np.tile(basis.local_dofs_index, len(basis.local_dofs_index))),
          axis=1)
      global_indices = np.stack(
          (np.repeat(dofs, len(dofs)), np.tile(dofs, len(dofs))), axis=1)
      # local_indices = get_indeces(basis.local_dofs_index)
      # global_indices = get_indeces(dofs)
      row = global_indices[:, 0]
      col = global_indices[:, 1]
      mat = self.elem_mat[i]
      if mat.TYPE in ['Fluid'] or (mat.TYPE in ['Poroelastic'] and 'P' in var):
        mat_coeff = 1 / (self.omega**2 * mat.rho_f)
      elif mat.TYPE in ['Poroelastic']:
        if var in ['Ux', 'Uy', 'Uz']:
          mat_coeff = mat.P_hat
      else:
        print("Material type not supported")
      elem_data = mat_coeff * basis.ke[local_indices[:, 0], local_indices[:,
                                                                          1]]
      rows.extend(row)
      cols.extend(col)
      data.extend(elem_data)
    self.K = coo_matrix((data, (rows, cols)),
                        shape=(self.nb_global_dofs, self.nb_global_dofs),
                        dtype=self.dtype).tocsr()

    return self.K

  def assemble_material_M(self, bases, var=None):
    if var is None:
      dofs_index = self.fe_space.get_global_dofs()
    else:
      dofs_index = self.fe_space.get_global_dofs_by_base(var)
    rows = []
    cols = []
    data = []
    for i, (dofs, basis) in enumerate(zip(dofs_index, bases)):
      local_indices = get_indeces(basis.local_dofs_index)
      global_indices = get_indeces(dofs)
      row = global_indices[:, 0]
      col = global_indices[:, 1]
      mat = self.elem_mat[i]
      if mat.TYPE in ['Fluid'] or (mat.TYPE in ['Poroelastic'] and 'P' in var):
        mat_coeff = 1 / mat.K_f
      elif mat.TYPE in ['Poroelastic']:
        if var in ['Ux', 'Uy', 'Uz']:
          mat_coeff = (self.omega**2) * mat.rho_til
      else:
        print("Material type not supported")
      elem_data = mat_coeff * basis.me[local_indices[:, 0], local_indices[:,
                                                                          1]]
      rows.extend(row)
      cols.extend(col)
      data.extend(elem_data)
    self.M = coo_matrix((data, (rows, cols)),
                        shape=(self.nb_global_dofs, self.nb_global_dofs),
                        dtype=self.dtype).tocsr()

    return self.M

  def assemble_global_material_matrix(self, bases, var=None):
    if var is None:
      dofs_index = self.fe_space.get_global_dofs()
    else:
      dofs_index = self.fe_space.get_global_dofs_by_base(var)
    rows = []
    cols = []
    data_K = []
    data_M = []
    for i, (dofs, basis) in enumerate(zip(dofs_index, bases)):
      local_indices = np.stack(
          (np.repeat(basis.local_dofs_index, len(basis.local_dofs_index)),
           np.tile(basis.local_dofs_index, len(basis.local_dofs_index))),
          axis=1)
      global_indices = np.stack(
          (np.repeat(dofs, len(dofs)), np.tile(dofs, len(dofs))), axis=1)
      row = global_indices[:, 0]
      col = global_indices[:, 1]
      mat = self.elem_mat[i]
      if mat.TYPE in ['Fluid'] or (mat.TYPE in ['Poroelastic'] and 'P' in var):
        mat_coeff_M = 1 / mat.K_f
        mat_coeff_K = 1 / (self.omega**2 * mat.rho_f)
      elif mat.TYPE in ['Poroelastic']:
        if var in ['Ux', 'Uy', 'Uz']:
          mat_coeff_M = (self.omega**2) * mat.rho_til
          mat_coeff_K = mat.P_hat
      else:
        print("Material type not supported")
      elem_data_M = mat_coeff_M * basis.me[local_indices[:, 0],
                                           local_indices[:, 1]]
      elem_data_K = mat_coeff_K * basis.ke[local_indices[:, 0],
                                           local_indices[:, 1]]
      rows.extend(row)
      cols.extend(col)
      data_K.extend(elem_data_K)
      data_M.extend(elem_data_M)
    self.M = coo_matrix((data_M, (rows, cols)),
                        shape=(self.nb_global_dofs, self.nb_global_dofs),
                        dtype=self.dtype).tocsr()
    self.K = coo_matrix((data_K, (rows, cols)),
                        shape=(self.nb_global_dofs, self.nb_global_dofs),
                        dtype=self.dtype).tocsr()

# ===================================== parallel assembly ==========================
# work not good

  def assembly_woker(self, args):
    dofs, basis, mat, omega, var, get_indeces = args
    local_indices = get_indeces(basis.local_dofs_index)
    global_indices = get_indeces(dofs)
    row = global_indices[:, 0]
    col = global_indices[:, 1]
    if mat.TYPE in ['Fluid'] or (mat.TYPE in ['Poroelastic'] and 'P' in var):
      mat_coeff_M = 1 / mat.K_f
      mat_coeff_K = 1 / (omega**2 * mat.rho_f)
    elif mat.TYPE in ['Poroelastic']:
      if var in ['Ux', 'Uy', 'Uz']:
        mat_coeff_M = (omega**2) * mat.rho_til
        mat_coeff_K = mat.P_hat
    else:
      print("Material type not supported")

    elem_data_M = mat_coeff_M * basis.me[local_indices[:, 0], local_indices[:,
                                                                            1]]
    elem_data_K = mat_coeff_K * basis.ke[local_indices[:, 0], local_indices[:,
                                                                            1]]
    return row, col, elem_data_M, elem_data_K

  def parallel_assemble_global_material_matrix(self, bases, var=None):
    if var is None:
      dofs_index = self.fe_space.get_global_dofs()
    else:
      dofs_index = self.fe_space.get_global_dofs_by_base(var)

    args = [(dofs, bases[i], self.elem_mat[i], self.omega, var, get_indeces)
            for i, dofs in enumerate(dofs_index)]
    # import pdb

    # for arg in args:
    # elems_data = self.assembly_woker(arg)
    # pdb.set_trace()
    with Pool() as pool:
      elems_data = pool.map(self.assembly_woker, args)

    rows = []
    cols = []
    data_M = []
    data_K = []

    for elem_data in elems_data:
      row, col, elem_data_M, elem_data_K = elem_data
      rows.extend(row)
      cols.extend(col)
      data_M.extend(elem_data_M)
      data_K.extend(elem_data_K)

    self.M = coo_matrix((data_M, (rows, cols)),
                        shape=(self.nb_global_dofs, self.nb_global_dofs),
                        dtype=self.dtype).tocsr()
    self.K = coo_matrix((data_K, (rows, cols)),
                        shape=(self.nb_global_dofs, self.nb_global_dofs),
                        dtype=self.dtype).tocsr()


# ===================================== end parallel assembly ==========================
class HelmholtzAssembler(BaseAssembler):

  def __init__(self, fe_space, omega, dtype) -> None:
    """
        General assembler for Helmholtz equation
        bases: list of basis
        subdomains: dict of subdomains
        dtype: data type of linear system"""
    super().__init__(fe_space, omega, dtype)
    self.elem_mat = {
        elem: mat
        for elem, mat in self.elem_mat.items() if mat.TYPE == 'Fluid'
    }
    [elem.set_frequency(omega) for elem in self.elem_mat.values()]

    self.K = csr_array((self.nb_global_dofs, self.nb_global_dofs),
                       dtype=self.dtype)
    self.M = csr_array((self.nb_global_dofs, self.nb_global_dofs),
                       dtype=self.dtype)

  def assembly_global_matrix(self, bases, var=None):
    self.initial_matrix()
    # self.assemble_material_K(bases, var)
    # self.assemble_material_M(bases, var)
    self.assemble_global_material_matrix(bases, var)
    # self.parallel_assemble_global_material_matrix(bases, var)

  def get_global_matrix(self):
    return self.K - self.M


class BiotAssembler(BaseAssembler):
  """
    Assembler for Biot's equation (only for Biot UP coupling equations)
    """

  def __init__(self, fe_space, omega, dtype) -> None:
    super().__init__(fe_space, omega, dtype)
    self.C = csr_array((self.nb_global_dofs, self.nb_global_dofs),
                       dtype=self.dtype)
    self.elem_mat = {
        elem: mat
        for elem, mat in self.elem_mat.items() if mat.TYPE == 'Poroelastic'
    }
    # update the material properties
    [elem.set_frequency(omega) for elem in self.elem_mat.values()]

    # self.elem_mat.update({i: mat for i in np.arange(len(elems))})

  def initial_matrix(self):
    self.K = csr_array((self.nb_global_dofs, self.nb_global_dofs),
                       dtype=self.dtype)
    self.M = csr_array((self.nb_global_dofs, self.nb_global_dofs),
                       dtype=self.dtype)
    self.C = csr_array((self.nb_global_dofs, self.nb_global_dofs),
                       dtype=self.dtype)

  def assemble_material_C(self, bases, var_1=None, var_2=None):
    if var_1 is None:
      dofs_index_1 = self.fe_space.get_global_dofs()
    else:
      dofs_index_1 = self.fe_space.get_global_dofs_by_base(var_1)
      dofs_index_2 = self.fe_space.get_global_dofs_by_base(var_2)
    for i, (dofs_1, dofs_2,
            basis) in enumerate(zip(dofs_index_1, dofs_index_2, bases)):
      local_indices = get_indeces(basis.local_dofs_index)
      global_indices = get_indeces(dofs_1, dofs_2)
      row = global_indices[:, 0]
      col = global_indices[:, 1]
      mat = self.elem_mat[i]
      if mat.TYPE in ['Poroelastic']:
        mat_coeff = mat.gamma_til
      else:
        print("Material type not supported")
      data = mat_coeff * basis.ce[local_indices[:, 0], local_indices[:, 1]]
      self.C += csr_array((data, (row, col)),
                          shape=(self.nb_global_dofs, self.nb_global_dofs),
                          dtype=self.dtype)

    return self.C

  def assembly_global_matrix(
      self,
      bases,
      vars,
  ):
    if len(bases) != len(vars) != 2:
      raise ValueError("the number of bases and variables have to be two")
    self.initial_matrix()
    K_p = self.assemble_material_K(bases[0], vars[0])
    M_p = self.assemble_material_M(bases[0], vars[0])
    self.initial_matrix()
    K_u = self.assemble_material_K(bases[1], vars[1])
    M_u = self.assemble_material_M(bases[1], vars[1])
    C_pu = self.assemble_material_C(bases[0], vars[1], vars[0])
    C_up = C_pu.T

    self.Mglobal = K_p + K_u - M_u - M_p - C_pu - C_up

  def get_global_matrix(self):
    return self.Mglobal


class CouplingAssember:
  """
    Assembly (combine) the (global) matrices of each component (Helmholz, elastic and Biot, etc)
    """

  def __init__(self,
               mesh,
               subdomains,
               components,
               coupling_type="continue") -> None:
    self.mesh = mesh
    self.subdomains = subdomains
    self.components = components
    self.coupling_type = coupling_type
    self.nb_global_dofs = 0
    self.nb_external_dofs = 0
    self.dtype = np.int8
    for comp in components:
      self.nb_global_dofs += comp.nb_global_dofs
      self.nb_external_dofs += comp.dof_handler.nb_external_dofs
      self.dtype = comp.dtype

    if "continue" in coupling_type:
      self.nb_global_dofs -= 1    # minus the duplicated continueous dofs
      self.nb_external_dofs -= 1
      self.nb_internal_dofs = self.nb_global_dofs - self.nb_external_dofs

    self.global_matrix = lil_array((self.nb_global_dofs, self.nb_global_dofs),
                                   dtype=self.dtype)

  def assembly_gloabl_matrix(self):
    # import pdb;pdb.set_trace()
    # first assembly the external dofs
    index_external_start = 0
    index_internal_start = self.nb_external_dofs
    for comp in self.components:
      local_external_index = comp.dof_handler.nb_external_dofs
      index_external_end = index_external_start + comp.dof_handler.nb_external_dofs
      self.global_matrix[
          index_external_start:index_external_end,
          index_external_start:index_external_end] += comp.get_global_matrix(
          ).tolil()[:local_external_index, :local_external_index]

      local_internal_index = local_external_index + comp.dof_handler.nb_internal_dofs
      index_internal_end = index_internal_start + comp.dof_handler.nb_internal_dofs
      self.global_matrix[
          index_internal_start:index_internal_end,
          index_internal_start:index_internal_end] += comp.get_global_matrix(
          ).tolil()[local_external_index:local_internal_index,
                    local_external_index:local_internal_index]

      index_external_start = index_external_end - 1
      index_internal_start = index_internal_end

    return self.global_matrix.tocsr()
