import numpy as np

from scipy.sparse import csr_array, csr_matrix, coo_matrix
from .Basis import Lobbato1DElement
from .Quadratures import GaussLegendreQuadrature
from .Polynomial import Lobatto

from .Polynomial import Lagrange2DTri
from .Quadratures import GaussLegendre2DTri


class ApplyBoundaryConditions:

  def __init__(self, mesh, FE_space, left_hand_side, right_hand_side, omega=0):
    self.mesh = mesh
    self.FE_space = FE_space
    self.elem_mat = FE_space.elem_mat
    self.left_hand_side = left_hand_side
    self.right_hand_side = right_hand_side
    self.nb_dofs = FE_space.get_nb_dofs()
    self.omega = omega
    self.dtype = self.right_hand_side.dtype

  def mesh2dof(self, position, var=None):
    return self.FE_space.get_dofs_from_var_coord(position, var)

  def apply_essential_bc(self,
                         essential_bcs,
                         var=None,
                         bctype='strong',
                         penalty=1e5):
    # import pdb; pdb.set_trace()
    dof_index = self.mesh2dof(essential_bcs['position'], var)

    if bctype == 'strong':
      left_hand_side_lil = self.left_hand_side.tolil()
      left_hand_side_lil[dof_index, :] = 0
      left_hand_side_lil[:, dof_index] = 0
      left_hand_side_lil[dof_index, dof_index] = 1
      self.left_hand_side = left_hand_side_lil.tocsr()
      self.right_hand_side[dof_index] += essential_bcs['value']

    elif bctype == 'penalty':
      row = np.array([dof_index])
      col = np.array([dof_index])
      mat = self.elem_mat[0]
      data = penalty * mat.P_hat * np.ones((1), dtype=self.dtype)
      # data = penalty*basis.me[1,1]
      self.left_hand_side += csr_array((data, (row, col)),
                                       shape=(self.nb_dofs, self.nb_dofs),
                                       dtype=self.dtype)
      self.right_hand_side[dof_index] += penalty * essential_bcs['value']

    elif bctype == 'nitsche':
      mat = self.elem_mat[0]
      alpha = 1e1
      scaling = 2 * self.mesh.get_nb_elems()
      nitsch = -1 * mat.P_hat * scaling * np.array(
          [[0, 0], [-0.5, 0.5]]) - 1 * mat.P_hat * scaling * np.array(
              [[0, -0.5], [0, 0.5]]) + alpha * np.array([[0, 0], [0, 1]])
      left_hand_side_lil = self.left_hand_side.tolil()
      left_hand_side_lil[dof_index - 1:dof_index + 1,
                         dof_index - 1:dof_index + 1] += nitsch
      self.left_hand_side = left_hand_side_lil.tocsr()

      self.right_hand_side[dof_index - 1] += 0.5 * essential_bcs['value']
      self.right_hand_side[dof_index] += alpha * essential_bcs[
          'value'] - 0.5 * essential_bcs['value']

    else:
      print("Weak imposing methods has not been implemented")

  def apply_source(self, source, bases, var=None):
    elements2node = self.mesh.get_mesh()
    lag2d_poly_o1 = Lagrange2DTri(1)
    points_o1, weights_o1 = GaussLegendre2DTri(3).points(), GaussLegendre2DTri(
        3).weights()
    if var is None:
      dofs_index = self.FE_space.get_global_dofs()
    else:
      dofs_index = self.FE_space.get_global_dofs_by_base(var)
    max_entries = len(dofs_index) * len(dofs_index[0]) * len(dofs_index[0])
    rows = np.empty(max_entries, dtype=int)
    data_M = np.empty(max_entries, dtype=self.dtype)

    idx = 0
    for i, (dofs, basis) in enumerate(zip(dofs_index, bases)):
      local_indices = basis.local_dofs_index
      global_indices = dofs

      node_coords = elements2node[i]
      xx = np.dot(basis.N, node_coords)
      f = np.array([
          weights_o1[i] * source['value'](x[0], x[1]) for i, x in enumerate(xx)
      ])
      integrand = basis.det_J * np.dot(basis.N.T,
                                       f) * 1 / (self.omega**2 * 1.213)

      elem_data_M = integrand[local_indices]
      size = len(global_indices)
      rows[idx:idx + size] = global_indices[:]
      data_M[idx:idx + size] = elem_data_M
      idx += size

    source_interp = csr_array((data_M[:idx], (rows[:idx], [0] * idx)),
                              shape=(self.nb_dofs, 1),
                              dtype=self.dtype)
    self.right_hand_side += source_interp
    return self.right_hand_side

  def apply_impedance_bc(self, impedence_bcs, var=None):
    if isinstance(impedence_bcs['position'], float):    #1D case
      dof_index = self.mesh2dof(impedence_bcs['position'], var)
      row = np.array([dof_index])
      col = np.array([dof_index])
      mat = self.elem_mat[dof_index - 1]
      mat.set_frequency(self.omega)
      mat_coeff = 1j * 1 / mat.rho_f * (self.omega /
                                        mat.c_f) * impedence_bcs['value']
      data = np.array([mat_coeff * 1])
      C_damp = csr_array((data, (row, col)),
                         shape=(self.nb_dofs, self.nb_dofs),
                         dtype=self.dtype)
      self.left_hand_side += C_damp
      return C_damp
    elif isinstance(impedence_bcs['position'], np.ndarray):
      edges = impedence_bcs['position']
      lines = self.mesh.exterior_facets[edges]
      rows = []
      cols = []
      data_M = []
      for line_index in lines:
        node_1_coord = self.mesh.nodes[line_index[0]]
        node_2_coord = self.mesh.nodes[line_index[1]]
        gl_q = GaussLegendreQuadrature(3)
        jac = np.linalg.norm(node_1_coord - node_2_coord) / 2
        gl_pts, gl_wts = gl_q.points(), gl_q.weights()
        l = Lobatto(1)
        N = l.get_shape_functions()
        f = np.zeros((2, 2), dtype=self.left_hand_side.dtype)
        for i, gl_pt in enumerate(gl_pts):
          x = N[0](gl_pt) * node_1_coord + N[1](gl_pt) * node_2_coord
          f_e = 1j / (self.omega * impedence_bcs['value'](x[0], x[1]))
          f += gl_wts[i] * np.array([[
              N[0](gl_pt) * N[0](gl_pt), N[0](gl_pt) * N[1](gl_pt)
          ], [N[1](gl_pt) * N[0](gl_pt), N[1](gl_pt) * N[1](gl_pt)]]) * f_e
        f *= jac    #
        local_index = np.arange(len(N))
        local_indices = np.stack(
            (np.repeat(local_index, len(N)), np.tile(local_index, len(N))),
            axis=1)
        global_indices = np.stack(
            (np.repeat(line_index,
                       len(line_index)), np.tile(line_index, len(line_index))),
            axis=1)
        elem_data_M = f[local_indices[:, 0], local_indices[:, 1]]
        row = global_indices[:, 0]
        col = global_indices[:, 1]
        rows.extend(row)
        cols.extend(col)
        data_M.extend(elem_data_M)

      C_damp = coo_matrix((data_M, (rows, cols)),
                          shape=(self.nb_dofs, self.nb_dofs),
                          dtype=self.dtype).tocsr()

      self.left_hand_side += C_damp
    return self.left_hand_side

  def apply_nature_bc(self, nature_bc, var=None, integr_order=1):
    if isinstance(nature_bc['position'], float):    #1D case
      dof_index = self.mesh2dof(nature_bc['position'], var)
      if nature_bc['type'] == 'fluid_velocity':
        self.right_hand_side[dof_index] += -nature_bc['value'] / (1j *
                                                                  self.omega)
      elif nature_bc['type'] == 'total_displacement':
        self.right_hand_side[dof_index] += nature_bc['value']
      elif nature_bc['type'] == 'solid_stress':
        self.right_hand_side[dof_index] += nature_bc['value']
      else:
        print("Nature BC type not supported")
    elif isinstance(nature_bc['position'], np.ndarray):
      breakpoint()
      edges = nature_bc['position']
      lines = self.mesh.exterior_facets[edges]
      for line_index in lines:
        compute_normal_vector(self.mesh, line_index)
        # breakpoint()
        gl_q = GaussLegendreQuadrature(10)
        gl_pts, gl_wts = gl_q.points(), gl_q.weights()
        l = Lobatto(1)
        N = l.get_shape_functions()
        basis = Lobbato1DElement(var, order=1, nodes=line_index)
        f = basis.integrate(lambda x: nature_bc['value'](x[0], x[1]), gl_pts,
                            gl_wts)
        # map coordiante into reference space
        if nature_bc['type'] == 'fluid_velocity':
          self.right_hand_side[line_index] += f / (1j * self.omega)
        elif nature_bc['type'] == 'analytical_gradient':
          self.right_hand_side[line_index] += f / self.omega**2
        elif nature_bc['type'] == 'total_displacement':
          self.right_hand_side[line_index] += f
        elif nature_bc['type'] == 'solid_stress':
          self.right_hand_side[line_index] += f
        else:
          print("Nature BC type not supported")


def compute_normal_vector(mesh, edge_or_facet):
  if len(edge_or_facet) == 2:
    edge = edge_or_facet
    node_1_coord = mesh.nodes[edge[0]]    # in physical space
    node_2_coord = mesh.nodes[edge[1]]
    tangent = np.array(
        [node_2_coord[0] - node_1_coord[0], node_2_coord[1] - node_1_coord[1]])
    normal = np.array([tangent[1], -tangent[0]])
    normal = normal / np.linalg.norm(normal)
  elif len(edge_or_facet) == 3:
    facet = edge_or_facet
    node_1_coord = mesh.nodes[facet[0]]    # in physical space
    node_2_coord = mesh.nodes[facet[1]]
    node_3_coord = mesh.nodes[facet[2]]
    vect_1 = node_2_coord - node_1_coord
    vect_2 = node_3_coord - node_1_coord

    # compute the normal vector with cross product
    normal = np.cross(vect_1, vect_2)
    normal = normal / np.linalg.norm(normal)
  return normal
