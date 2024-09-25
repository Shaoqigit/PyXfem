import os

current_dir = os.path.dirname(os.path.realpath(__file__))
working_dir = os.path.join(current_dir, "..")
import sys

sys.path.append(working_dir)

from SAcouS.acxfem import Lagrange3DTetraElement, Lagrange2DTriElement
import numpy as np


def construct_3D_tetra(order):
  label = "fluid"
  if order == 1:
    nodes = np.array([[0, 0, 0], [1 * 2, 0, 0], [0, 1 * 3, 0], [0, 0, 1 / 2]])
  lag_3d_tetra = Lagrange3DTetraElement(label, order, nodes)
  return lag_3d_tetra


def construct_2D_tri(order):
  label = "fluid"
  if order == 2:
    nodes = np.array([[0, 0], [1, 0], [0, 1], [1 / 2, 0], [1 / 2, 1 / 2],
                      [0, 1 / 2]])
  lag_2d_tri = Lagrange2DTriElement(label, order, nodes)
  return lag_2d_tri


def reference_matrices(dim, order):
  if dim == 3:
    if order == 1:
      Ke_ref = np.array([[0.5, -0.16666667, -0.16666667, -0.16666667],
                         [-0.16666667, 0.16666667, 0., 0.],
                         [-0.16666667, 0., 0.16666667, 0.],
                         [-0.16666667, 0., 0., 0.16666667]])
      Me_ref = np.array([[0.016666667, 0.008333333, 0.008333333, 0.008333333],
                         [0.008333333, 0.016666667, 0.008333333, 0.008333333],
                         [0.008333333, 0.008333333, 0.016666667, 0.008333333],
                         [0.008333333, 0.008333333, 0.008333333, 0.016666667]])
    elif order == 2:
      raise NotImplementedError
  elif dim == 2:
    if order == 2:
      Ke_ref = np.array(
          [[1.0, 0.166666667, 0.166666667, -0.666666667, -0.666666667, 0.],
           [0.166666667, 0.5, 0., -0.666666667, 0., 0.],
           [0.166666667, 0., 0.5, 0., -0.666666667, 0.],
           [-0.666666667, -0.666666667, 0., 2.666666667, 0., -1.33333333],
           [-0.666666667, 0., -0.666666667, 0., 2.666666667, -1.33333333],
           [0., 0., 0., -1.33333333, -1.333333333, 2.666666667]])
      Me_ref = np.array([[
          1.66666667e-02, -2.77777778e-03, -2.77777778e-03,
          6.76542156e-17 - 6.24500451e-17, -1.11111111e-02
      ],
                         [
                             -2.77777778e-03, 1.66666667e-02, -2.77777778e-03,
                             -3.94649591e-17, -1.11111111e-02, 7.11507676e-17
                         ],
                         [
                             -2.77777778e-03, -2.77777778e-03, 1.66666667e-02,
                             -1.11111111e-02, 9.36750677e-17, -5.67850888e-17
                         ],
                         [
                             6.76542156e-17, -3.94649591e-17, -1.11111111e-02,
                             8.88888889e-02, 4.44444444e-02, 4.44444444e-02
                         ],
                         [
                             -6.24500451e-17, -1.11111111e-02, 9.36750677e-17,
                             4.44444444e-02, 8.88888889e-02, 4.44444444e-02
                         ],
                         [
                             -1.11111111e-02, 7.11507676e-17, -5.67850888e-17,
                             4.44444444e-02, 4.44444444e-02, 8.88888889e-02
                         ]])
    elif order == 1:
      raise NotImplementedError
  return Ke_ref, Me_ref


def test_case():
  # plot the nodes
  lag_3d_tetra = construct_3D_tetra(1)
  Ke_ref, Me_ref = reference_matrices(3, 1)

  lag_2d_tri_p2 = construct_2D_tri(2)
  print(lag_2d_tri_p2.ke)
  print(lag_2d_tri_p2.me)
  # print(lag_3d_tetra.ke)
  # print(lag_3d_tetra.me)
  #compare the results
  if np.allclose(lag_3d_tetra.ke, Ke_ref) and np.allclose(
      lag_3d_tetra.me, Me_ref):
    print("Test passed!")
    return True
  else:
    print("Test failed!")
    return False


if __name__ == "__main__":
  result = test_case()
