from skfem import *
from skfem.helpers import dot, grad
from skfem.visuals.matplotlib import draw
from SAcouS.acxfem.mesh import Mesh2D
import numpy as np

mesh = Mesh2D()
mesh.read_mesh("mesh/square_1.msh")
# # enable additional mesh validity checks, sacrificing performance
# import logging
# logging.basicConfig(format='%(levelname)s %(asctime)s %(name)s %(message)s')
# logging.getLogger('skfem').setLevel(logging.DEBUG)

# create the mesh
# read mesh from gmsh .mesh file
# read the mesh from a .mshe file
# or, with your own points and cells:
points = mesh.nodes
cells = mesh.connectivity
m = MeshTri(points.T, cells.T)
# plot the mesh
draw(m)
e = ElementTriP1()
basis = Basis(m, e)
facebasis = FacetBasis(m, e, facets=m.boundaries['left'])

freq = 2000
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
  return 1 * np.exp(-1j * omega)


@LinearForm
def numann_flux(v, w):
  x, y = w.x
  return neumann_bc(x, y) * v


A = asm(laplace, basis)
M = asm(mass, basis)
b = asm(numann_flux, facebasis) * omega * 1.213
# or:
# A = laplace.assemble(basis)
# b = rhs.assemble(basis)
# import pdb

# pdb.set_trace()

# # pdb.set_trace()
# print(A.toarray())
# print(M.toarray())
K = A - k**2 * M
from scipy.sparse.linalg import spsolve

u = spsolve(K, b)
# solve -- can be anything that takes a sparse matrix and a right-hand side
x = solve(K, b)


def visualize():
  from skfem.visuals.matplotlib import plot
  return plot(m, x, shading='gouraud', colorbar=True)


if __name__ == "__main__":
  visualize().show()
