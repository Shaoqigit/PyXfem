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

# utilities.py
# including some useful functions that are independent of the other modules
from scipy.sparse import csr_array

def check_material_compability(subdomains):
    mats = []
    for key in subdomains.keys():
        mats.append(key)

    compatibale_mats = mats[0].COMPATIBLE
    print(compatibale_mats)
    for mat in mats:
        print(mat.TYPE)
        if mat.TYPE not in compatibale_mats:
            raise ValueError("Material model is not compatible")
        else:
            print("Material models are compatible, computation continues ...")


def display_matrix_in_array(M:csr_array):
        return M.toarray()


def plot_matrix_partten(M:csr_array):
    import matplotlib.pyplot as plt
    plt.spy(M)
    plt.show()