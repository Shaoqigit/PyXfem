import numpy as np
from src.quadratures import GaussLegendreQuadrature
from src.polynomial import Lobatto, Larange
from numpy.polynomial.legendre import leggauss

def compute_matrix(Ke, Me, order):
    n_pts = order*2
    gl_q = GaussLegendreQuadrature(n_pts)
    gl_pts, gl_wts = gl_q.points(), gl_q.weights()
    gl_pts, gl_wts = leggauss(n_pts)
    l = Lobatto(order)
    B = l.get_der_shape_functions()
    N = l.get_shape_functions()
    len(Me[0])
    for i in range(len(Me[0])):
        for j in range(len(Me[0])):
            Ke[i, j] = sum(gl_wt*B[i](gl_pt)*B[j](gl_pt) for gl_pt, gl_wt in zip(gl_pts, gl_wts))
            Me[i, j] = sum(gl_wt*N[i](gl_pt)*N[j](gl_pt) for gl_pt, gl_wt in zip(gl_pts, gl_wts))
            if abs(Ke[i, j]) < 1e-14:
                Ke[i, j] = 0
            if abs(Me[i, j]) < 1e-14:
                Me[i, j] = 0
            


# 1D lobatto element matrix: p=1
order = 1
Ke1Do1 = np.zeros((2,2))
Me1Do1 = np.zeros((2,2))
compute_matrix(Ke1Do1, Me1Do1, order)


# 1D lobatto element matrix: p=2
order = 2
Ke1Do2 = np.zeros((3,3))
Me1Do2 = np.zeros((3,3))
compute_matrix(Ke1Do2, Me1Do2, order)

# 1D lobatto element  matrix: p=3
order = 3
Ke1Do3 = np.zeros((4,4))
Me1Do3 = np.zeros((4,4))
compute_matrix(Ke1Do3, Me1Do3, order)

# 1D lobatto element matrix: p=4
order = 4
Ke1Do4 = np.zeros((5,5))
Me1Do4 = np.zeros((5,5))
compute_matrix(Ke1Do4, Me1Do4, order)

Ke1D = [Ke1Do1, Ke1Do2, Ke1Do3, Ke1Do4]
Me1D = [Me1Do1, Me1Do2, Me1Do3, Me1Do4]

# print(Ke1Do1)
# print(Me1Do1)
# print(Ke1Do2)
# print(Me1Do2)
# print(Ke1Do3)
# print(Me1Do3)
# print(Ke1Do4)
# print(Me1Do4)