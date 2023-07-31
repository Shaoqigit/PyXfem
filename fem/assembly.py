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
    
    def assemble_material_C(self, omega = 0):
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
    
    
    def assemble_nature_bc(self, nature_bc):
        F = np.zeros(self.num_dofs, dtype=self.dtype)
        if nature_bc['type']=='velocity':
            F[nature_bc['position']] = 1j * self.omega * nature_bc['value']
        else:
            print("Nature BC type not supported")

        return F
    
    def get_dim(self):
        return self.num_dofs


class Assembler4Biot:
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

        

    def assemble_material_K(self, bases, var = None, omega = 0):
        self.omega = omega
        if var is None:
            dofs_index = self.dof_handler.get_global_dofs()
        else:
            dofs_index = self.dof_handler.get_global_dofs_by_base(var)

        for i, (dofs, basis) in enumerate(zip(dofs_index, bases)):
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
            elif mat.TYPE in ['Poroelastic']:
                mat.set_frequency(omega)
                if var == 'P':
                    mat_coeff = 1/(self.omega**2*mat.rho_f)
                elif var in ['Ux', 'Uy', 'Uz']:
                    mat_coeff = 1/mat.P_hat
            else:
                print("Material type not supported")
            data = mat_coeff*basis.ke[local_indices.T[0], local_indices.T[1]]
            self.K += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs), dtype=self.dtype)

        return self.K
    
    def assemble_material_M(self, bases, var = None, omega = 0):
        self.omega = omega
        if var is None:
            dofs_index = self.dof_handler.get_global_dofs()
        else:
            dofs_index = self.dof_handler.get_global_dofs_by_base(var)
        for i, (dofs, basis) in enumerate(zip(dofs_index, bases)):
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
            elif mat.TYPE in ['Poroelastic']:
                mat.set_frequency(omega)
                if var == 'P':
                    mat_coeff = 1/mat.rho_f*(1/mat.c_f)**2
                elif var in ['Ux', 'Uy', 'Uz']:
                    mat_coeff = (omega**2)*mat.rho_s_til
            else:
                print("Material type not supported")
            data = mat_coeff*basis.me[local_indices.T[0], local_indices.T[1]]
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
            local_indices = np.array([(row, col) for row in basis.local_dofs_index() for col in basis.local_dofs_index()])
            global_indices = np.array([(row, col) for row in dofs_1 for col in dofs_2])
            # print(global_indices)
            row = global_indices.T[0]
            col = global_indices.T[1]
            mat = self.elem_mat[i]
            if mat.TYPE in ['Air', 'Fluid']:
                mat_coeff = 1/mat.rho_f*(omega/mat.c_f)**2
            elif mat.TYPE in ['Equivalent Fluid', 'Limp Fluid']:
                mat.set_frequency(omega)
                mat_coeff = 1/mat.rho_f*(omega/mat.c_f)**2
            elif mat.TYPE in ['Poroelastic']:
                mat.set_frequency(omega)
                mat_coeff = mat.gamma_til
            else:
                print("Material type not supported")
            data = mat_coeff*basis.ce[local_indices.T[0], local_indices.T[1]]
            self.C += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs), dtype=self.dtype)

        return self.C
    
    def apply_essential_bc(self, left_hand_side, essential_bcs, var=None, type='strong'):
        if type == 'strong':
            left_hand_side.tolil()
            left_hand_side[essential_bcs['position'], :] = 0
            left_hand_side[:, essential_bcs['position']] = 0
            left_hand_side[essential_bcs['position'], essential_bcs['position']] = 1
            left_hand_side.tocsr()
        else:
            print("Weak imposing methods has not been implemented")

    def apply_impedance_bc(self, impedence_bcs):
        row = np.array([impedence_bcs['position']])
        col = np.array([impedence_bcs['position']])
        mat = self.elem_mat[impedence_bcs['position']-1]
        mat.set_frequency(self.omega)
        mat_coeff = 1j*1/mat.rho_f*(self.omega/mat.c_f)*impedence_bcs['value']
        data = np.array([mat_coeff*1])
        C = csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs), dtype=self.dtype)

        return C
    
    
    def apply_nature_bc(self, nature_bc):
        F = np.zeros(self.num_dofs, dtype=self.dtype)
        if nature_bc['type']=='velocity':
            F[nature_bc['position']] = 1j * self.omega * nature_bc['value']
        elif nature_bc['type']=='total_displacement':
            F[nature_bc['position']] = nature_bc['value']
        else:
            print("Nature BC type not supported")

        return F
    
    def get_dim(self):
        return self.num_dofs
