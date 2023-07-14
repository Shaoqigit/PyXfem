import numpy as np
from abc import ABCMeta, abstractmethod
from scipy.sparse import csc_matrix

from scipy.sparse.linalg import spsolve, bicg, bicgstab, gmres, SuperLU

from src.polynomial import Lobatto, Larange
from src.quadratures import GaussLegendreQuadrature

class BaseSolver(metaclass=ABCMeta):
    """base abstract FE solver class

    """
    def __init__(self, dof_handler):
        self.internal_dofs = dof_handler.num_internal_dofs
        self.external_dofs = dof_handler.num_external_dofs
        self.u = None

    @abstractmethod
    def solve(self):
        pass

class LinearSolver(BaseSolver):
    """linear solver class
    parameters:
    left_hand_side: ndarray
        left hand side matrix
    right_hand_side: ndarray
        right hand side vector
    """
    def solve(self, left_hand_side, right_hand_side):
        u = spsolve(left_hand_side, right_hand_side)
        self.u = u[:self.external_dofs]

    def condensation(self):
        """condensation of linear system
        parameters:
        internal_dofs: ndarray
            internal dofs
        """

        K_ee = self.lhs[:self.external_dofs, :self.external_dofs]
        K_ei = self.lhs[:self.external_dofs, self.external_dofs:]
        K_ei = self.lhs[:self.external_dofs, self.external_dofs:]
        K_ii = self.lhs[self.external_dofs:, self.external_dofs:]

        schur_complement = K_ee - K_ei.dot(np.linalg.inv(K_ii)).dot(K_ei.T)
        self.lhs = self.lhs[self.internal_dofs, :][:, self.internal_dofs]
        self.rhs = self.rhs[self.internal_dofs]

    def simple_condensation
    