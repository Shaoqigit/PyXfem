import numpy as np
from numba import jit
from scipy.sparse import csr_array, lil_array

from acxtmm.adm_basis import AdmFluid as fluid_elem
from acxtmm.adm_basis import AdmPoroElastic2 as poroelastic_elem

class AdmAssembler:
    def __init__(self, mesh, subdomains, omega, dtype):
        self.dtype = dtype
        self.omega = omega
        self.mesh = mesh.get_mesh()
        self.nb_dofs = mesh.get_nb_nodes()
        self.elem_mats = {}
        for key, elems in subdomains.items():
            self.elem_mats.update({elem: key for elem in elems})


    def assemble_global_adm(self, theta, k_0, mode='continue'):
        """assemble admittance matrix
        parameters:
        mat: admittance material
        """
        self.global_adm = lil_array((self.nb_dofs, self.nb_dofs), dtype=self.dtype)
        for i, elem in self.mesh.items():
            mat = self.elem_mats[i]
            mat.set_frequency(self.omega)
            if mat.TYPE in ['Fluid']:
                adm = fluid_elem(mat, self.omega, theta, k_0, elem, mode)
                adm.admittance()
            elif mat.TYPE in ['Poroelastic']:
                adm = poroelastic_elem(mat, self.omega, theta, k_0, elem, mode)
                adm.admittance()
            self.global_adm[i:i+2, i:i+2] += adm.adm
            
        return self.global_adm.tocsr()
    
    def assemble_nature_bc(self, nature_bc):
        F = np.zeros(self.nb_dofs, dtype=self.dtype)
        if nature_bc['type']=='velocity':
            F[nature_bc['position']] = -1*nature_bc['value']/(1j*self.omega)
        elif nature_bc['type'] == 'total_displacement':
            F[nature_bc['position']] = nature_bc['value']
        else:
            print("Nature BC type not supported")

        return F