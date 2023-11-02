import numpy as np
from acxfem.solver import BaseSolver
from scipy import sparse
from scipy.sparse.linalg import spsolve


class EigenSolver(BaseSolver):
    """eigen solver class
    parameters:
    left_hand_side: ndarray
        left hand side matrix
    right_hand_side: ndarray
        right hand side vector
    """
    def solve(self, stiffness_matrix, mass_matrix, nb_modes):
        import time
        start = time.time()
        from scipy.linalg import eigh

        end = time.time()
        print("Eigenvalue problem solving time: ", end-start)

        # eig_freq, modes = eigs(stiffness_matrix, k=nb_modes, M=mass_matrix, which='LM', return_eigenvectors=True)
        eig_freq, modes = eigh(stiffness_matrix.toarray(), mass_matrix.toarray(), type=1, subset_by_index=[0, nb_modes-1])

        return eig_freq, modes
    
    
class ModalReduction:
    """modal reduction class
    parameters:
    K_w: ndarray original stiffness matrix
    M_w: ndarray original mass matrix
    modes: ndarray modal matrix [solution at modes]
    """
    def __init__(self, K_w, M_w, modes):
        self.K_w = K_w
        self.M_w = M_w
        self.Phi_m = modes

    def projection(self, original_array):
        """project original matrix/vector on modal basis to modal space
        parameters:
        original_array: ndarray
        returns:
        array: ndarray
            projected reducded matrix/vector
        """
        if len(original_array.shape) == 1:
            array = self.Phi_m.T @ original_array
            return array
        else:
            array = self.Phi_m.T @ original_array @ self.Phi_m
            return sparse.csr_matrix(array)
    
    def solve(self, left_hand_matrix, right_hand_vector):
        """solve the reduced system
        parameters:
        left_hand_matrix: ndarray
            left hand matrix
        right_hand_vector: ndarray
            right hand vector
        returns:
        sol: ndarray
            solution
        """
        import time
        start = time.time()
        sol = spsolve(left_hand_matrix, right_hand_vector)
        end = time.time()
        print("Reduced Linear solving time: ", end-start)
        return sol
    
    def recover_sol(self, reducde_sol):
        """recover the physical solution from modal space
        parameters:
        reduced_sol: ndarray/solution vector

        returns:
        sol: recoverd solution in physical space
        """
        sol = self.Phi_m @ reducde_sol
        return sol
    