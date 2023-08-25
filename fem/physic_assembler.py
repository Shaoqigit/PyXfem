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

# assembly the global/partial matrices according to the physic of the components
from fem.mesh import Mesh1D
from fem.dofhandler import DofHandler1D

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

class BaseAssembler:
    def __init__(self, dof_handler, subdomains, dtype) -> None:
        """
        General assembler for Helmholtz equation
        bases: list of basis
        subdomains: dict of subdomains
        dtype: data type of linear system"""
        self.dof_handler = dof_handler
        self.num_dofs = dof_handler.get_num_dofs()
        self.dtype = dtype
        self.elem_mat = {}
        for mat, elems in subdomains.items():
            self.elem_mat.update({i: mat for i in np.arange(len(elems))})

        self.K = csr_array((self.num_dofs, self.num_dofs), dtype=self.dtype)
        self.M = csr_array((self.num_dofs, self.num_dofs), dtype=self.dtype)

    def initial_matrix(self):
        self.K = csr_array((self.num_dofs, self.num_dofs), dtype=self.dtype)
        self.M = csr_array((self.num_dofs, self.num_dofs), dtype=self.dtype)

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
                if 'P' in var:
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
                if 'P' in var:
                    mat_coeff = 1/mat.K_eq_til
                elif var in ['Ux', 'Uy', 'Uz']:
                    mat_coeff = (omega**2)*mat.rho_til
            else:
                print("Material type not supported")
            data = mat_coeff*basis.me[local_indices[:,0], local_indices[:,1]]
            self.M += csr_array((data, (row, col)), shape=(self.num_dofs, self.num_dofs), dtype=self.dtype)

        return self.M


class HelmholtzAssembler(BaseAssembler):
    def __init__(self, dof_handler, subdomains, dtype) -> None:
        """
        General assembler for Helmholtz equation
        bases: list of basis
        subdomains: dict of subdomains
        dtype: data type of linear system"""
        super().__init__(dof_handler, subdomains, dtype)
        self.elem_mat = {}
        for mat, elems in subdomains.items():
            if mat.TYPE == 'Fluid':
                self.elem_mat.update({elem: mat for elem in elems})

        self.K = csr_array((self.num_dofs, self.num_dofs), dtype=self.dtype)
        self.M = csr_array((self.num_dofs, self.num_dofs), dtype=self.dtype)
        

class BiotAssembler(BaseAssembler):
    """
    Assembler for Biot's equation (only for Biot UP coupling equations)
    """
    def __init__(self, dof_handler, subdomains, dtype) -> None:
        super().__init__(dof_handler, subdomains, dtype)
        import pdb; pdb.set_trace()
        self.C = csr_array((self.num_dofs, self.num_dofs), dtype=self.dtype)
        self.elem_mat = {}
        for mat, elems in subdomains.items():
            if mat.TYPE == 'Poroelastic':
                self.elem_mat.update({i: mat for i in np.arange(len(elems))})


    def initial_matrix(self):
        self.K = csr_array((self.num_dofs, self.num_dofs), dtype=self.dtype)
        self.M = csr_array((self.num_dofs, self.num_dofs), dtype=self.dtype)
        self.C = csr_array((self.num_dofs, self.num_dofs), dtype=self.dtype)


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
    
# class CouplingAssember: