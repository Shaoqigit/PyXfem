import numpy as np
from fem.solver import BaseSolver


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
        # from scipy.linalg import eig
        # from scipy.sparse.linalg import eigsh
        # from scipy.sparse.linalg import eigs
        from scipy.linalg import eigh

        end = time.time()
        print("Eigenvalue problem solving time: ", end-start)

        # eig_freq, modes = eigs(stiffness_matrix, k=nb_modes, M=mass_matrix, which='LM', return_eigenvectors=True)
        eig_freq, modes = eigh(stiffness_matrix.toarray(), mass_matrix.toarray(), type=1, subset_by_index=[0, nb_modes-1])

        return eig_freq, modes
    
    
class ModalReduction:
    """modal reduction class
    parameters:
    order: int
        element order
    nodes: ndarray
        1d: [x1, x2]
        2d: [(x1, y1), (x2, y2)]
        3d: [(x1, y1, z1), (x2, y2, z2)]
    """
    def __init__(self, eigenvalue, modes):
        self.eigenvalues = eigenvalue
        self.Phi_m = modes

    def projection(self, original_array):
        """project f to modal space
        parameters:
        f: ndarray
            f(x)
        returns:
        f_m: ndarray
            f_m(x)
        """
        if len(original_array.shape) == 1:
            array = self.Phi_m.T @ original_array
        else:
            array = self.Phi_m.T @ original_array @ self.Phi_m
        return array
    
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
        sol = np.linalg.solve(left_hand_matrix, right_hand_vector)
        return sol
    
    def recover_sol(self, reducde_sol):
        """recover the solution from modal space
        parameters:
        f_m: ndarray
            f_m(x)
        returns:
        f: ndarray
            f(x)
        """
        sol = self.Phi_m @ reducde_sol
        return sol
    