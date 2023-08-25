import numpy as np
from numba import jit
from scipy.sparse import csr_array, lil_array


@jit(nopython=True)
def get_indeces(*dofs):
    if len(dofs) == 1:    
        return np.array([(row, col) for row in dofs[0] for col in dofs[0]])
    elif len(dofs) == 2:
        return np.array([(row, col) for row in dofs[0] for col in dofs[1]])
    else:
        raise ValueError("wrong number of arguments")

# @jit(nopython=True)
# def general_K_assembler(matrix, num_dofs, dtype, dofs_index, bases, elem_mats, omega, var = None):
#     for i, (dofs, basis) in enumerate(zip(dofs_index, bases)):
#         local_indices = get_indeces(basis.local_dofs_index())
#         global_indices = get_indeces(dofs)
#         # print(global_indices)
#         row = global_indices[:,0]
#         col = global_indices[:,1]
#         mat = elem_mats[i]
#         if mat.TYPE in ['Air', 'Fluid']:
#             mat_coeff = 1/mat.rho_f
#         elif mat.TYPE in ['Equivalent Fluid', 'Limp Fluid']:
#             mat.set_frequency(omega)
#             mat_coeff = 1/mat.rho_f
#         elif mat.TYPE in ['Poroelastic']:
#             mat.set_frequency(omega)
#             if var == 'P':
#                 mat_coeff = 1/(omega**2*mat.rho_f)
#             elif var in ['Ux', 'Uy', 'Uz']:
#                 mat_coeff = mat.P_hat
#         else:
#             print("Material type not supported")
#         data = mat_coeff*basis.ke[local_indices[:,0], local_indices[:,1]]
#         matrix[row, col] += data
#     return matrix

class Assembler:
    def __init__(self, dof_handler, bases, subdomains, dtype) -> None:
        """
        dof_handler: DofHandler1D
        bases: list of basis
        subdomains: dict of subdomains
        dtype: data type of linear system"""

        self.dof_handler = dof_handler
        self.num_dofs = dof_handler.get_num_dofs()
        self.bases = bases
        self.dtype = dtype
        self.K = csr_array((self.num_dofs, self.num_dofs), dtype=self.dtype)
        self.M = csr_array((self.num_dofs, self.num_dofs), dtype=self.dtype)
        self.omega = 0.
        self.elem_mat = {}
        for key, elems in subdomains.items():
            self.elem_mat.update({elem: key for elem in elems})

    def assemble_K(self):
        """
        get the global stiffness matrix without material property
        """

        for dofs, basis in zip(self.dof_handler.get_global_dofs(), self.bases):
            local_indices = np.array([(row, col) for row in basis.local_dofs_index() for col in basis.local_dofs_index()])
            global_indices = np.array([(row, col) for row in dofs for col in dofs])
            # print(global_indices)
            row = global_indices.T[0]
            col = global_indices.T[1]
            data = basis.ke[local_indices.T[0], local_indices.T[1]]
            self.K += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs))

    def assemble_M(self): 
        """
        get the global mass matrix without material property
        """

        for dofs, basis in zip(self.dof_handler.get_global_dofs(), self.bases):
            local_indices = np.array([(row, col) for row in basis.local_dofs_index() for col in basis.local_dofs_index()])
            global_indices = np.array([(row, col) for row in dofs for col in dofs])
            # print(global_indices)
            row = global_indices.T[0]
            col = global_indices.T[1]
            data = basis.me[local_indices.T[0], local_indices.T[1]]
            self.M += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs))

        

    def assemble_material_K(self, omega = 0):
        """
        get the global stiffness matrix with material property (frequency dependent material)
        """
        self.omega = omega
        for i, (dofs, basis) in enumerate(zip(self.dof_handler.get_global_dofs(), self.bases)):
            local_indices = get_indeces(basis.local_dofs_index())
            global_indices = get_indeces(dofs)
            # print(global_indices)
            row = global_indices[:,0]
            col = global_indices[:,1]
            mat = self.elem_mat[i]
            mat.set_frequency(omega)
            if mat.TYPE in ['Fluid']:
                mat_coeff = 1/mat.rho_f
            else:
                print("Material type not supported")
            data = mat_coeff*basis.ke[local_indices.T[0], local_indices[:,1]]
            self.K += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs), dtype=self.dtype)

        return self.K
    
    def assemble_material_M(self, omega = 0):
        """
        get the global mass matrix with material property (frequency dependent material)
        """
        self.omega = omega
        for i, (dofs, basis) in enumerate(zip(self.dof_handler.get_global_dofs(), self.bases)):
            local_indices = get_indeces(basis.local_dofs_index())
            global_indices = get_indeces(dofs)
            # print(global_indices)
            row = global_indices[:,0]
            col = global_indices[:,1]
            mat = self.elem_mat[i]
            mat.set_frequency(omega)
            if mat.TYPE in ['Fluid']:
                mat_coeff = 1/mat.rho_f*(omega/mat.c_f)**2
            else:
                print("Material type not supported")
            data = mat_coeff*basis.me[local_indices[:,0], local_indices[:,1]]
            self.M += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs), dtype=self.dtype)

        return self.M
    
    def assemble_material_C(self, omega = 0):
        """
        get the global coupling matrix with material property (frequency dependent material)
        """
        self.omega = omega
        for i, (dofs, basis) in enumerate(zip(self.dof_handler.get_global_dofs(), self.bases)):
            local_indices = np.array([(row, col) for row in basis.local_dofs_index() for col in basis.local_dofs_index()])
            global_indices = np.array([(row, col) for row in dofs for col in dofs])
            # print(global_indices)
            row = global_indices[:,0]
            col = global_indices[:,1]
            mat = self.elem_mat[i]
            mat.set_frequency(omega)
            if mat.TYPE in ['Fluid']:
                mat_coeff = 1/mat.rho_f*(omega/mat.c_f)**2
            else:
                print("Material type not supported")
            data = mat_coeff*basis.me[local_indices[:,0], local_indices[:,1]]
            self.M += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs), dtype=self.dtype)

        return self.M
    

    def assemble_impedance_bc(self, impedence_bcs):
        """
        impedence_bcs: dict of impedence boundary conditions
        return the global matrix contributed from impedence boundary conditions
        """
        row = np.array([impedence_bcs['position']])
        col = np.array([impedence_bcs['position']])
        mat = self.elem_mat[impedence_bcs['position']-1]
        mat.set_frequency(self.omega)
        mat_coeff = 1j*1/mat.rho_f*(self.omega/mat.c_f)*impedence_bcs['value']
        data = np.array([mat_coeff*1])
        C = csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs), dtype=self.dtype)

        return C
    
    
    def assemble_nature_bc(self, nature_bc):
        """
        nature_bc: dict of nature boundary conditions
        return the global vector contributed from nature boundary conditions
        """
        F = np.zeros(self.num_dofs, dtype=self.dtype)
        if nature_bc['type']=='velocity':
            F[nature_bc['position']] = 1j * self.omega * nature_bc['value']
        else:
            print("Nature BC type not supported")

        return F
    
    def get_dim(self):
        return self.num_dofs


class Assembler4Biot:
    """
    Assembler for Biot's equation (only for Biot UP coupling equations)
    """
    def __init__(self, dof_handler, subdomains, dtype) -> None:
        self.dof_handler = dof_handler
        self.num_dofs = dof_handler.get_num_dofs()
        self.dtype = dtype
        self.K = csr_array((self.num_dofs, self.num_dofs), dtype=self.dtype)
        self.M = csr_array((self.num_dofs, self.num_dofs), dtype=self.dtype)
        self.C = csr_array((self.num_dofs, self.num_dofs), dtype=self.dtype)
        self.omega = 0.
        self.elem_mat = {}
        for key, elems in subdomains.items():
            self.elem_mat.update({elem: key for elem in elems})

    def assemble_K(self, bases):
        for dofs, basis in zip(self.dof_handler.get_global_dofs(), bases):
            local_indices = np.array([(row, col) for row in basis.local_dofs_index() for col in basis.local_dofs_index()])
            global_indices = np.array([(row, col) for row in dofs for col in dofs])
            # print(global_indices)
            row = global_indices.T[0]
            col = global_indices.T[1]
            data = basis.ke[local_indices.T[0], local_indices.T[1]]
            self.K += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs))

    def assemble_M(self, bases):  
        for dofs, basis in zip(self.dof_handler.get_global_dofs(), bases):
            local_indices = np.array([(row, col) for row in basis.local_dofs_index() for col in basis.local_dofs_index()])
            global_indices = np.array([(row, col) for row in dofs for col in dofs])
            # print(global_indices)
            row = global_indices.T[0]
            col = global_indices.T[1]
            data = basis.me[local_indices.T[0], local_indices.T[1]]
            self.M += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs))

    def initial_matrix(self):
        self.K = csr_array((self.num_dofs, self.num_dofs), dtype=self.dtype)
        self.M = csr_array((self.num_dofs, self.num_dofs), dtype=self.dtype)
        self.C = csr_array((self.num_dofs, self.num_dofs), dtype=self.dtype)

    def assemble_material_K(self, bases, var = None, omega = 0):
        self.omega = omega
        if var is None:
            dofs_index = self.dof_handler.get_global_dofs()
        else:
            dofs_index = self.dof_handler.get_global_dofs_by_base(var)
        for i, (dofs, basis) in enumerate(zip(dofs_index, bases)):
            local_indices = get_indeces(basis.local_dofs_index())
            global_indices = get_indeces(dofs)
            row = global_indices[:,0]
            col = global_indices[:,1]
            mat = self.elem_mat[i]
            mat.set_frequency(omega)
            if mat.TYPE in ['Fluid']:
                mat_coeff = 1/mat.rho_f
            elif mat.TYPE in ['Poroelastic']:
                if var == 'P':
                    mat_coeff = 1/(self.omega**2*mat.rho_f)
                elif var in ['Ux', 'Uy', 'Uz']:
                    mat_coeff = mat.P_hat
            else:
                print("Material type not supported")
            data = mat_coeff*basis.ke[local_indices[:,0], local_indices[:,1]]
            self.K += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs), dtype=self.dtype)

        return self.K


    
    def assemble_material_M(self, bases, var = None, omega = 0):
        self.omega = omega
        if var is None:
            dofs_index = self.dof_handler.get_global_dofs()
        else:
            dofs_index = self.dof_handler.get_global_dofs_by_base(var)
        for i, (dofs, basis) in enumerate(zip(dofs_index, bases)):
            local_indices = get_indeces(basis.local_dofs_index())
            global_indices = get_indeces(dofs)
            # print(global_indices)
            row = global_indices[:,0]
            col = global_indices[:,1]
            mat = self.elem_mat[i]
            mat.set_frequency(omega)
            if mat.TYPE in ['Fluid']:
                mat_coeff = 1/mat.rho_f*(omega/mat.c_f)**2
            elif mat.TYPE in ['Poroelastic']:
                mat.set_frequency(omega)
                if var == 'P':
                    mat_coeff = 1/mat.K_eq_til
                elif var in ['Ux', 'Uy', 'Uz']:
                    mat_coeff = (omega**2)*mat.rho_til
            else:
                print("Material type not supported")
            data = mat_coeff*basis.me[local_indices[:,0], local_indices[:,1]]
            self.M += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs), dtype=self.dtype)

        return self.M
    

    def assemble_material_C(self, bases, var_1=None, var_2=None, omega = 0):
        self.omega = omega
        if var_1 is None:
            dofs_index_1 = self.dof_handler.get_global_dofs()
        else:
            dofs_index_1 = self.dof_handler.get_global_dofs_by_base(var_1)
            dofs_index_2 = self.dof_handler.get_global_dofs_by_base(var_2)
        for i, (dofs_1, dofs_2, basis) in enumerate(zip(dofs_index_1, dofs_index_2, bases)):
            local_indices = get_indeces(basis.local_dofs_index())
            global_indices = get_indeces(dofs_1, dofs_2)
            # print(global_indices)
            row = global_indices[:,0]
            col = global_indices[:,1]
            mat = self.elem_mat[i]
            mat.set_frequency(omega)
            if mat.TYPE in ['Poroelastic']:
                mat_coeff = mat.gamma_til
            else:
                print("Material type not supported")
            data = mat_coeff*basis.ce[local_indices[:,0], local_indices[:,1]]
            self.C += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs), dtype=self.dtype)

        return self.C
    
    def apply_essential_bc(self, left_hand_side, essential_bcs, var=None, bctype='strong', penalty=1e5):
        dof_index = self.dof_handler.mesh2dof(essential_bcs['position'], var)

        if bctype == 'strong':
            left_hand_side_lil = left_hand_side.tolil()
            left_hand_side_lil[dof_index, :] = 0
            left_hand_side_lil[:, dof_index] = 0
            left_hand_side_lil[dof_index, dof_index] = 1
            left_hand_side = left_hand_side_lil.tocsr()
            # left_hand_side
            # left_hand_side[dof_index, :] = 0
            # left_hand_side[:, dof_index] = 0
            # left_hand_side[dof_index, dof_index] = 1
            self.F = None
            return left_hand_side

        elif bctype == 'penalty':
            row = np.array([dof_index])
            col = np.array([dof_index])
            mat = self.elem_mat[0]
            data = penalty*mat.P_hat*np.ones((1), dtype=self.dtype)
            # data = penalty*basis.me[1,1]
            left_hand_side += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs), dtype=self.dtype)
            self.F = np.zeros(self.num_dofs, dtype=self.dtype)
            self.F[dof_index] = penalty*essential_bcs['value']
            return left_hand_side
        elif bctype == 'nitsche':
            mat = self.elem_mat[0]
            alpha = 1e1
            scaling = 2*1000
            nitsch = -1*mat.P_hat*scaling*np.array([[0, 0], [-0.5, 0.5]])-1*mat.P_hat*scaling*np.array([[0, -0.5], [0, 0.5]])+alpha*np.array([[0, 0], [0, 1]])
            left_hand_side_lil = left_hand_side.tolil()
            left_hand_side_lil[dof_index-1:dof_index+1, dof_index-1:dof_index+1] += nitsch
            left_hand_side = left_hand_side_lil.tocsr()

            self.F = np.zeros(self.num_dofs, dtype=self.dtype)
            self.F[dof_index-1] += 0.5*essential_bcs['value']
            self.F[dof_index] += alpha*essential_bcs['value']-0.5*essential_bcs['value']
            return left_hand_side
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
        C = csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs), dtype=self.dtype)

        return C
    
    
    def apply_nature_bc(self, nature_bc, var=None):
        dof_index = self.dof_handler.mesh2dof(nature_bc['position'], var)
        if self.F is None:
            right_hand_side = np.zeros(self.num_dofs, dtype=self.dtype)
        else:
            right_hand_side = self.F
        if nature_bc['type']=='velocity':
            right_hand_side[dof_index] += 1j * self.omega * nature_bc['value']
        elif nature_bc['type']=='total_displacement':
            right_hand_side[dof_index] += nature_bc['value']
        else:
            print("Nature BC type not supported")

        return right_hand_side
    
    def get_dim(self):
        return self.num_dofs
