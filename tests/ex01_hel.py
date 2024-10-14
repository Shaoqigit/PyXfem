from skfem import *
from skfem.helpers import dot, grad
from skfem.visuals.matplotlib import draw
from SAcouS.Mesh import Mesh2D, MeshReader

import numpy as np
import os
# mesh.read_mesh("mesh/square_1.msh")
# mesh_reader = MeshReader("mesh/square_1.msh")
# current_dir = os.path.dirname(os.path.realpath(__file__))
# mesh_reader = MeshReader(current_dir + "/mesh/unit_tube_3D_refine.msh", dim=3)
# mesh = mesh_reader.get_mesh()
# # enable additional mesh validity checks, sacrificing performance
# import logging
# logging.basicConfig(format='%(levelname)s %(asctime)s %(name)s %(message)s')
# logging.getLogger('skfem').setLevel(logging.DEBUG)

# create the mesh
# read mesh from gmsh .mesh file
# read the mesh from a .mshe file
# or, with your own points and cells:
# points = mesh.nodes
# cells = mesh.elem_connect
# points = np.array([[0, 0, 0], [1 * 2, 0, 0], [0, 1 * 3, 0], [0, 0, 1 / 2]])
# cells = np.array([[0, 1, 2, 3]])
points = np.array([[0, 0], [1, 0], [0, 1], [1 / 2, 0], [1 / 2, 1 / 2],
                   [0, 1 / 2]])
cells = np.array([[0, 1, 2, 3, 4, 5]])
# point
# m = MeshTet(points.T, cells.T)
m = MeshTri2(points.T, cells.T)
# plot the mesh
# draw(m)
e = ElementTriP2()
# e = ElementTetP1()
basis = Basis(m, e)

# facebasis = FacetBasis(m, e)

freq = 200
omega = 2 * np.pi * freq
k = omega / 343


# this method could also be imported from skfem.models.laplace
@BilinearForm
def laplace(u, v, _):
  return dot(grad(u), grad(v))


@BilinearForm
def mass(u, v, _):
  return u * v


# this method could also be imported from skfem.models.unit_load
def neumann_bc(x, y):
  return 1


@LinearForm
def numann_flux(v, w):
  x, y = w.x
  return neumann_bc(x, y) * v


# A = asm(laplace, basis) / (omega**2 * 1.213)
K = asm(laplace, basis)
# M = asm(mass, basis) / (1.400 * 1.01325e5)
M = asm(mass, basis)
print(K.toarray())
print(M.toarray())
# b = asm(numann_flux, facebasis) / (omega * 1j) * np.exp(-1j * omega)

import pdb

pdb.set_trace()
# breakpoint()
# or:
# A = laplace.assemble(basis)
# b = rhs.assemble(basis)
# import pdb

# pdb.set_trace()

# # pdb.set_trace()
# print(A.toarray())
# print(M.toarray())
# K = A - M
# from scipy.sparse.linalg import spsolve

# u = spsolve(K, b)
# # solve -- can be anything that takes a sparse matrix and a right-hand side
# x = solve(K, b)


def visualize():
  from skfem.visuals.matplotlib import plot
  return plot(m, x.real, shading='gouraud', colorbar=True)


if __name__ == "__main__":
  visualize().show()
