import os

current_dir = os.path.dirname(os.path.realpath(__file__))
working_dir = os.path.join(current_dir, "..")
import sys

sys.path.append(working_dir)

from SAcouS.acxfem import Lagrange3DTetraElement
import numpy as np


def test_case():
  label = "fluid"
  order = 1
  nodes = np.array([[0, 0, 0], [1 * 2, 0, 0], [0, 1 * 3, 0], [0, 0, 1 / 2]])
  # plot the nodes
  lag_3d_tetra = Lagrange3DTetraElement(label, order, nodes)
  import matplotlib.pyplot as plt
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  ax.scatter(nodes[:, 0], nodes[:, 1], nodes[:, 2])
  plt.show()

  Ke_ref = np.array([[0.5, -0.16666667, -0.16666667, -0.16666667],
                     [-0.16666667, 0.16666667, 0., 0.],
                     [-0.16666667, 0., 0.16666667, 0.],
                     [-0.16666667, 0., 0., 0.16666667]])
  Me_ref = np.array([[0.016666667, 0.008333333, 0.008333333, 0.008333333],
                     [0.008333333, 0.016666667, 0.008333333, 0.008333333],
                     [0.008333333, 0.008333333, 0.016666667, 0.008333333],
                     [0.008333333, 0.008333333, 0.008333333, 0.016666667]])
  print(lag_3d_tetra.ke)
  print(lag_3d_tetra.me)
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
