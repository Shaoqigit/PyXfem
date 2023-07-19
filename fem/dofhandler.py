# This file is part of PyXfem, a software distributed under the MIT license.
# For any question, please contact the authors cited below.
#
# Copyright (c) 2023
# 	Shaoqi WU <shaoqiwu@outlook.com
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

# dofhandle.py
# deal with mesh and basis, return global dof index t

import numpy as np
from fem.mesh import Mesh1D



class DofHandler1D:

    def __init__(self, mesh, bases) -> None:
        if isinstance(mesh, Mesh1D):
            self.mesh = mesh
            self.bases = bases
            self.connect = mesh.connectivity
        else:
            raise TypeError("mesh must be 1D mesh")


    def get_num_dofs(self):
        """return global dof"""
        num_dofs = 0
        for basis in self.bases:
            if basis.is_discontinue:
                num_discontiue = self.basis.interface
                return self.mesh.get_num_elems()*self.basis.get_order() +1 + num_discontiue*(self.basis.get_order()+1)
            else:
                num_dofs += basis.get_order()
        return num_dofs+1

            
    
    def get_global_dofs(self):
        """return local dof
        return: global dof index
        [node_1_index, node_2_index, internal_dof_index],
        [.....],
        [.....]]"""
        global_dof = []
        internal_dof_index_start = self.mesh.get_num_nodes()
        for i, basis in enumerate(self.bases):
            # print("internal_dof_index_start", internal_dof_index_start)
            if basis.num_internal_dofs() == 0:
                global_dof.append(self.mesh.connectivity[i])
                continue
            elem = self.mesh.connectivity[i]
            global_dof.append(np.hstack((np.array([elem[0], elem[1]]), np.array([internal_dof_index_start + j for j in range(0, basis.get_order()-1)]))))
            internal_dof_index_start += basis.get_order()-1

        return global_dof
        
    @property
    def num_external_dofs(self):
        return self.mesh.get_num_nodes()
    
    @property
    def num_internal_dofs(self):
        return self.get_num_dofs() - self.num_external_dofs

