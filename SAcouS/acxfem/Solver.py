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

# solver is used solve and optimize the linear system

import numpy as np
from abc import ABCMeta, abstractmethod
from scipy.sparse import csr_array

from scipy.sparse.linalg import spsolve, bicg, bicgstab, gmres, SuperLU
from scipy.sparse.csgraph import reverse_cuthill_mckee

from .Polynomial import Lobatto, Larange

try:
  from petsc4py import PETSc
  from petsc4py.PETSc import Mat, Vec
  PETSC_on = True
except ImportError:
  PETSC_on = False
  pass


class BaseSolver(metaclass=ABCMeta):
  """base abstract FE solver class

    """

  def __init__(self, fe_space=None, coupling_assember=None, symmetric=True):
    if fe_space is None:
      self.internal_dofs = coupling_assember.nb_internal_dofs
      self.external_dofs = coupling_assember.nb_external_dofs
      self.nb_dofs = coupling_assember.nb_global_dofs
    else:
      self.internal_dofs = fe_space.nb_internal_dofs
      self.external_dofs = fe_space.nb_external_dofs
      self.nb_dofs = fe_space.nb_dofs
    if fe_space is None and coupling_assember is None:
      raise ValueError(
          "dof_handler and coupling_assember cannot be None at the same time")

    self.u = None
    self.sym = symmetric

  @abstractmethod
  def solve(self):
    pass


def petsc_solver(left_hand_side, right_hand_side):
  """petsc solver
    parameters:
    left_hand_side: CSR matrix
        left hand side matrix
    right_hand_side: ndarray
        right hand side vector
    """
  if not isinstance(left_hand_side, PETSc.Mat):
    left_hand_side = PETSc.Mat().createAIJ(size=left_hand_side.shape,
                                           csr=(left_hand_side.indptr,
                                                left_hand_side.indices,
                                                left_hand_side.data))
    left_hand_side.assemblyBegin()
    left_hand_side.assemblyEnd()
  else:
    assert (isinstance(left_hand_side, PETSc.Mat))
  b = PETSc.Vec().createWithArray(right_hand_side)
  x = PETSc.Vec().createSeq(right_hand_side.shape[0])
  ksp = PETSc.KSP().create()
  ksp.setOperators(left_hand_side)
  ksp.setType(
      PETSc.KSP.Type.PREONLY)    # Use preconditioner only (direct solve)

  pc = ksp.getPC()
  pc.setType(PETSc.PC.Type.LU)    # Set the preconditioner to LU
  pc.setFactorSolverType(PETSc.Mat.SolverType.MUMPS)

  ksp.setFromOptions()
  ksp.solve(b, x)
  print('PETSc dirct solver used')
  return x.getArray()


class LinearSolver(BaseSolver):
  """linear solver class
    parameters:
    left_hand_side: ndarray
        left hand side matrix
    right_hand_side: ndarray
        right hand side vector
    """

  def solve(self, left_hand_side, right_hand_side, solver='spsolve'):
    import time
    start = time.time()
    if solver == 'petsc' and PETSC_on:
      u = petsc_solver(left_hand_side, right_hand_side)
    else:
      u = spsolve(left_hand_side, right_hand_side)
    # u = np.linalg.solve(left_hand_side.toarray(), right_hand_side)
    self.u = u[:self.external_dofs]
    end = time.time()
    print("Linear system solving time: ", end - start)

  def condition_number(self, left_hand_side):
    return np.linalg.cond(left_hand_side.toarray())

  def optimize_matrix_pattern(self, left_hand_side, right_hand_side):
    row = left_hand_side.tocoo().row
    col = left_hand_side.tocoo().col
    perm_rcm = reverse_cuthill_mckee(left_hand_side, symmetric_mode=True)
    iperm_rcm = np.zeros(shape=self.nb_dofs, dtype=np.int32)
    for js, jt in enumerate(perm_rcm):
      iperm_rcm[jt] = js
    I_rcm = iperm_rcm[row]
    J_rcm = iperm_rcm[col]
    V_rcm = left_hand_side.data
    print(I_rcm)
    left_hand_side = csr_array((V_rcm, (I_rcm, J_rcm)),
                               shape=(self.nb_dofs, self.nb_dofs))
    # import pdb; pdb.set_trace()
    right_hand_side[I_rcm[0]], right_hand_side[row[0]] = right_hand_side[
        row[0]], right_hand_side[I_rcm[0]]
    return left_hand_side, right_hand_side

  def static_condensation(self):
    """condensation of linear system
        parameters:
        internal_dofs: ndarray
            internal dofs
        """

    K_ee = self.lhs[:self.external_dofs, :self.external_dofs]
    K_ei = self.lhs[:self.external_dofs, self.external_dofs:]
    K_ie = self.lhs[self.external_dofs:, :self.external_dofs]
    K_ii = self.lhs[self.external_dofs:, self.external_dofs:]

    schur_complement = K_ee - K_ei.dot(np.linalg.inv(K_ii)).dot(K_ei.T)
    self.rhs = K_ei.dot(np.linalg.inv(K_ii)).dot(self.rhs[self.external_dofs:])
    return schur_complement, self.rhs


class AdmittanceSolver:
  """admittance solver class
    parameters:
    admittance: ndarray
        admittance matrix
    right_hand_side: ndarray
        right hand side vector
    """

  def __init__(self, admittance, right_hand_side):
    self.admittance = admittance
    self.right_hand_side = right_hand_side
    self.sol = None

  def solve(self):
    u = spsolve(self.admittance, self.right_hand_side)
    self.sol = u
