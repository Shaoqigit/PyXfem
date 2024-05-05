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

# Main Test case 
import os
current_dir = os.path.dirname(os.path.realpath(__file__))
working_dir = os.path.join(current_dir , "..")
import sys
sys.path.append(working_dir)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import spy

from SAcouS.acxfem.basis import Lobbato1DElement
from SAcouS.acxfem.mesh import Mesh1D
import meshio
from SAcouS.acxfem.dofhandler import DofHandler1D
from SAcouS.acxfem.assembly import Assembler
from SAcouS.acxfem.materials import Air, Fluid, EquivalentFluid
from SAcouS.acxfem.utilities import check_material_compability, display_matrix_in_array, plot_matrix_partten
from SAcouS.acxfem.solver import LinearSolver
from SAcouS.acxfem.postprocess import PostProcessField
from analytical.fluid_sol import ImpedenceKundltTube


def test_case_2D():
    mesh = meshio.read('mesh/tube_c.msh')
    connectivity = mesh.cells_dict['triangle']
    nodes = mesh.points
    from SAcouS.acxfem.precompute_matrices_lag import Ke1D, Me1D, Ce1D, compute_matrix





if __name__ == "__main__":
    result = test_case_2D()