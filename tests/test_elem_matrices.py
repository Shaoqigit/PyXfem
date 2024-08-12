import os

current_dir = os.path.dirname(os.path.realpath(__file__))
working_dir = os.path.join(current_dir, "..")
import sys

sys.path.append(working_dir)

from SAcouS.acxfem import Lagrange2DTriElement, Lagrange3DTetraElement
import numpy as np


def test_case():
  label = "fluid"
  order = 1

  correct_matrices = []

  # 2D p1 triangular element
  nodes = np.array([[0, 0], [2, -1], [1, 0.5]])
  lag_2d_tri = Lagrange2DTriElement(label, order, nodes)
  Ke_ref = np.array([[0.8125, 0.0625, -0.875], [0.0625, 0.3125, -0.375],
                     [-0.875, -0.375, 1.25]])
  Me_ref = np.array([[0.16666667, 0.08333333, 0.08333333],
                     [0.08333333, 0.16666667, 0.08333333],
                     [0.08333333, 0.08333333, 0.16666667]])
  if np.allclose(lag_2d_tri.ke, Ke_ref) and np.allclose(lag_2d_tri.me, Me_ref):
    correct_matrices.append(True)
  #compare the
  nodes = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])

  # 3D p1 tetrahedral element
  lag_3d_tetra = Lagrange3DTetraElement(label, order, nodes)

  Ke_ref = np.array([[0.5, -0.16666667, -0.16666667, -0.16666667],
                     [-0.16666667, 0.16666667, 0., 0.],
                     [-0.16666667, 0., 0.16666667, 0.],
                     [-0.16666667, 0., 0., 0.16666667]])
  Me_ref = np.array([[0.016666667, 0.008333333, 0.008333333, 0.008333333],
                     [0.008333333, 0.016666667, 0.008333333, 0.008333333],
                     [0.008333333, 0.008333333, 0.016666667, 0.008333333],
                     [0.008333333, 0.008333333, 0.008333333, 0.016666667]])

  if np.allclose(lag_3d_tetra.ke, Ke_ref) and np.allclose(
      lag_3d_tetra.me, Me_ref):
    correct_matrices.append(True)
  if all(correct_matrices):
    print("Test passed!")
    return True
  else:
    print("Test failed!")
    return False


if __name__ == "__main__":
  result = test_case()
