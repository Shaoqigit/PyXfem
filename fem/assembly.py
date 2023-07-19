import numpy as np
from scipy.sparse import csr_array, csr_matrix


class Assembler:
    def __init__(self, dof_handler, bases, subdomains, dtype) -> None:
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
        for dofs, basis in zip(self.dof_handler.get_global_dofs(), self.bases):
            local_indices = np.array([(row, col) for row in basis.local_dofs_index() for col in basis.local_dofs_index()])
            global_indices = np.array([(row, col) for row in dofs for col in dofs])
            # print(global_indices)
            row = global_indices.T[0]
            col = global_indices.T[1]
            data = basis.ke[local_indices.T[0], local_indices.T[1]]
            self.K += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs))

    def assemble_M(self):  
        for dofs, basis in zip(self.dof_handler.get_global_dofs(), self.bases):
            local_indices = np.array([(row, col) for row in basis.local_dofs_index() for col in basis.local_dofs_index()])
            global_indices = np.array([(row, col) for row in dofs for col in dofs])
            # print(global_indices)
            row = global_indices.T[0]
            col = global_indices.T[1]
            data = basis.me[local_indices.T[0], local_indices.T[1]]
            self.M += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs))

        

    def assemble_material_K(self, omega = 0):
        self.omega = omega
        for i, (dofs, basis) in enumerate(zip(self.dof_handler.get_global_dofs(), self.bases)):
            local_indices = np.array([(row, col) for row in basis.local_dofs_index() for col in basis.local_dofs_index()])
            global_indices = np.array([(row, col) for row in dofs for col in dofs])
            # print(global_indices)
            row = global_indices.T[0]
            col = global_indices.T[1]
            mat = self.elem_mat[i]
            if mat.TYPE in ['Air', 'Fluid']:
                mat_coeff = 1/mat.rho_f
            elif mat.TYPE in ['Equivalent Fluid', 'Limp Fluid']:
                mat.set_frequency(omega)
                mat_coeff = 1/mat.rho_f
            else:
                print("Material type not supported")
            data = mat_coeff*basis.ke[local_indices.T[0], local_indices.T[1]]
            self.K += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs), dtype=self.dtype)

        return self.K
    
    def assemble_material_M(self, omega = 0):
        self.omega = omega

        for i, (dofs, basis) in enumerate(zip(self.dof_handler.get_global_dofs(), self.bases)):
            local_indices = np.array([(row, col) for row in basis.local_dofs_index() for col in basis.local_dofs_index()])
            global_indices = np.array([(row, col) for row in dofs for col in dofs])
            # print(global_indices)
            row = global_indices.T[0]
            col = global_indices.T[1]
            mat = self.elem_mat[i]
            if mat.TYPE in ['Air', 'Fluid']:
                mat_coeff = 1/mat.rho_f*(omega/mat.c_f)**2
            elif mat.TYPE in ['Equivalent Fluid', 'Limp Fluid']:
                mat.set_frequency(omega)
                mat_coeff = 1/mat.rho_f*(omega/mat.c_f)**2
            else:
                print("Material type not supported")
            data = mat_coeff*basis.me[local_indices.T[0], local_indices.T[1]]
            self.M += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs), dtype=self.dtype)

        return self.M
    

    def assemble_impedance_bc(self, impedence_bcs):
        row = np.array([impedence_bcs['position']])
        col = np.array([impedence_bcs['position']])
        mat = self.elem_mat[impedence_bcs['position']-1]
        mat.set_frequency(self.omega)
        mat_coeff = 1j*1/mat.rho_f*(self.omega/mat.c_f)*impedence_bcs['value']
        data = np.array([mat_coeff*1])
        C = csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs), dtype=self.dtype)

        return C
    
    def get_dim(self):
        return self.num_dofs
    
    
    def assemble_nature_bc(self, nature_bc):
        F = np.zeros(self.num_dofs, dtype=self.dtype)
        if nature_bc['type']=='velocity':
            F[nature_bc['position']] = 1j * self.omega * nature_bc['value']
            # F = csr_matrix(F, shape=(self.num_dofs, 1), dtype=self.dtype)
        else:
            print("Nature BC type not supported")

        return F

