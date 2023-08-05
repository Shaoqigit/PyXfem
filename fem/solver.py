import numpy as np
from abc import ABCMeta, abstractmethod
from scipy.sparse import csr_array

from scipy.sparse.linalg import spsolve, bicg, bicgstab, gmres, SuperLU
from scipy.sparse.csgraph import reverse_cuthill_mckee

from fem.polynomial import Lobatto, Larange
from fem.quadratures import GaussLegendreQuadrature

class BaseSolver(metaclass=ABCMeta):
    """base abstract FE solver class

    """
    def __init__(self, dof_handler):
        self.internal_dofs = dof_handler.num_internal_dofs
        self.external_dofs = dof_handler.num_external_dofs
        self.num_dofs = dof_handler.get_num_dofs()
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
        # u = np.linalg.solve(left_hand_side.toarray(), right_hand_side)
        self.u = u[:self.external_dofs]

    def condition_number(self, left_hand_side):
        return np.linalg.cond(left_hand_side.toarray())

    def optimize_matrix_pattern(self, left_hand_side, right_hand_side):
        row = left_hand_side.tocoo().row
        col = left_hand_side.tocoo().col
        perm_rcm = reverse_cuthill_mckee(left_hand_side,symmetric_mode=True)
        iperm_rcm = np.zeros(shape=self.num_dofs,dtype=np.int32)
        for js,jt in enumerate(perm_rcm):
            iperm_rcm[jt] = js
        I_rcm   = iperm_rcm[row]
        J_rcm   = iperm_rcm[col]
        V_rcm   = left_hand_side.data
        print(I_rcm)
        left_hand_side = csr_array((V_rcm,(I_rcm,J_rcm)),shape=(self.num_dofs,self.num_dofs))
        import pdb; pdb.set_trace()
        right_hand_side[I_rcm[0]],right_hand_side[row[0]] = right_hand_side[row[0]],right_hand_side[I_rcm[0]]
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
    