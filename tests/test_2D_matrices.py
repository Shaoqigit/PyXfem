# import os

# current_dir = os.path.dirname(os.path.realpath(__file__))
# working_dir = os.path.join(current_dir, "..")
# import sys

# sys.path.append(working_dir)

from SAcouS.acxfem.basis import Lagrange2DTriElement
import numpy as np


def test_case():
  label = "fluid"
  order = 1
  nodes = np.array([[0, 0], [2, -1], [1, 0.5]])
  # plot the nodes
  lag_2d_tri = Lagrange2DTriElement(1, order, nodes)
  Ke_ref = np.array([[0.8125, 0.0625, -0.875], [0.0625, 0.3125, -0.375],
                     [-0.875, -0.375, 1.25]])
  Me_ref = np.array([[0.16666667, 0.08333333, 0.08333333],
                     [0.08333333, 0.16666667, 0.08333333],
                     [0.08333333, 0.08333333, 0.16666667]])
  print(lag_2d_tri.ke)
  print(lag_2d_tri.me)
  #compare the results
  if np.allclose(lag_2d_tri.ke, Ke_ref) and np.allclose(lag_2d_tri.me, Me_ref):
    print("Test passed!")
    return True
  else:
    print("Test failed!")
    return False


if __name__ == "__main__":
  result = test_case()
