import numpy as np
from numba import jit
from scipy.sparse import csr_array, lil_array

class ApplyBoundaryConditions:
    def __init__(self, mesh, FE_space, left_hand_side, right_hand_side, omega=0):
        self.mesh = mesh
        self.FE_space = FE_space
        self.elem_mat = FE_space.elem_mat
        self.left_hand_side = left_hand_side
        self.right_hand_side = right_hand_side
        self.nb_dofs = self.left_hand_side.shape[0]
        self.omega = omega

    def mesh2dof(self, position, var=None):
        return self.FE_space.get_dofs_from_var_coord(position, var)

    def apply_essential_bc(self, essential_bcs, var=None, bctype='strong', penalty=1e5):
        # import pdb; pdb.set_trace()
        dof_index = self.mesh2dof(essential_bcs['position'], var)

        if bctype == 'strong':
            left_hand_side_lil = self.left_hand_side.tolil()
            left_hand_side_lil[dof_index, :] = 0
            left_hand_side_lil[:, dof_index] = 0
            left_hand_side_lil[dof_index, dof_index] = 1
            self.left_hand_side = left_hand_side_lil.tocsr()
            self.right_hand_side[dof_index] += essential_bcs['value']

        elif bctype == 'penalty':
            row = np.array([dof_index])
            col = np.array([dof_index])
            mat = self.elem_mat[0]
            data = penalty*mat.P_hat*np.ones((1), dtype=self.dtype)
            # data = penalty*basis.me[1,1]
            self.left_hand_side += csr_array((data, (row, col)), shape=(self.nb_dofs, self.nb_dofs), dtype=self.dtype)
            self.right_hand_side[dof_index] += penalty*essential_bcs['value']

        elif bctype == 'nitsche':
            mat = self.elem_mat[0]
            alpha = 1e1
            scaling = 2*self.mesh.get_nb_elems()
            nitsch = -1*mat.P_hat*scaling*np.array([[0, 0], [-0.5, 0.5]])-1*mat.P_hat*scaling*np.array([[0, -0.5], [0, 0.5]])+alpha*np.array([[0, 0], [0, 1]])
            left_hand_side_lil = self.left_hand_side.tolil()
            left_hand_side_lil[dof_index-1:dof_index+1, dof_index-1:dof_index+1] += nitsch
            self.left_hand_side = left_hand_side_lil.tocsr()

            self.right_hand_side[dof_index-1] += 0.5*essential_bcs['value']
            self.right_hand_side[dof_index] += alpha*essential_bcs['value']-0.5*essential_bcs['value']

        else: 
            print("Weak imposing methods has not been implemented")

    def apply_impedance_bc(self, impedence_bcs, var=None):
        dof_index = self.dof_handler.mesh2dof(impedence_bcs['position'], var)
        row = np.array([dof_index])
        col = np.array([dof_index])
        mat = self.elem_mat[dof_index-1]
        mat.set_frequency(self.omega)
        mat_coeff = 1j*1/mat.rho_f*(self.omega/mat.c_f)*impedence_bcs['value']
        data = np.array([mat_coeff*1])
        self.left_hand_side += csr_array((data, (row, col)), shape=(self.nb_dofs, self.nb_dofs), dtype=self.dtype)

        
    def apply_nature_bc(self, nature_bc, var=None):
        dof_index = self.mesh2dof(nature_bc['position'], var)

        if nature_bc['type']=='velocity':
            self.right_hand_side[dof_index] += -nature_bc['value']/(1j*self.omega)
        elif nature_bc['type']=='total_displacement':
            self.right_hand_side[dof_index] += nature_bc['value']
        elif nature_bc['type']=='solid_stress':
            self.right_hand_side[dof_index] += nature_bc['value']
        else:
            print("Nature BC type not supported")

