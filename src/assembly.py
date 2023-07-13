import numpy as np
from scipy.sparse import csr_array, csr_matrix, csc_matrix, coo_matrix


class Assembler:
    def __init__(self, dof_handler, bases) -> None:
        self.dof_handler = dof_handler
        self.num_dofs = dof_handler.get_num_dofs()
        self.bases = bases
        self.M = csr_array((self.num_dofs, self.num_dofs))

    def assemble_K(self):
        for dofs, basis in zip(self.dof_handler.global_dof_index(), self.bases):
            local_indices = np.array([(row, col) for row in basis.local_dofs_index() for col in basis.local_dofs_index()])
            global_indices = np.array([(row, col) for row in dofs for col in dofs])
            # print(global_indices)
            row = global_indices.T[0]
            col = global_indices.T[1]
            data = basis.ke[local_indices.T[0], local_indices.T[1]]
            self.M += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs))

    def assemble_M(self):  
        for dofs, basis in zip(self.dof_handler.global_dof_index(), self.bases):
            local_indices = np.array([(row, col) for row in basis.local_dofs_index() for col in basis.local_dofs_index()])
            global_indices = np.array([(row, col) for row in dofs for col in dofs])
            # print(global_indices)
            row = global_indices.T[0]
            col = global_indices.T[1]
            data = basis.me[local_indices.T[0], local_indices.T[1]]
            self.M += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs))

        

    def assemble_material_K(self, subdomains):
        elem_mat = {}
        for key, elems in subdomains.items():
            elem_mat.update({elem: key for elem in elems})

        for i, dofs, basis in enumerate(zip(self.dof_handler.global_dof_index(), self.bases)):
            local_indices = np.array([(row, col) for row in basis.local_dofs_index() for col in basis.local_dofs_index()])
            global_indices = np.array([(row, col) for row in dofs for col in dofs])
            # print(global_indices)
            row = global_indices.T[0]
            col = global_indices.T[1]

            data = basis.ke[local_indices.T[0], local_indices.T[1]]
            self.M += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs))

        return self.M
    
    def get_matrix_in_array(self):
        return self.M.toarray()

