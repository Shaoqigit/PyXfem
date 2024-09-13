from collections import defaultdict
import numpy as np

from SAcouS.interface.parser import ParserFactory

from SAcouS.Mesh import Mesh1D
from SAcouS.PostProcess import PostProcessField, PostProcessFRF

from SAcouS.acxfem import Lobbato1DElement
from SAcouS.acxfem import FESpace
from SAcouS.acxfem import HelmholtzAssembler
from SAcouS.acxfem import ApplyBoundaryConditions
from SAcouS.acxfem import LinearSolver
from SAcouS.acxfem import check_material_compability


class PyAcoustiXSetuper:

  def __init__(self):
    self.file_path = None
    self.solution_pool = None

  def welcome(self):
    print(
        f'Welcome to Simple Acoustic Simulator based on lib PyAcoustiX!\nAuthor: Shaoqiwu@outlook.com\n'
    )
    print(f'Parsing the input file: {self.file_path}')

  def parse_input(self, file_path):
    self.file_path = file_path
    parser = ParserFactory.create_parser(self.file_path)
    self.sol_info = parser.parse()

  def exit(self):
    print('now exiting PyAcoustiX...')

  def fem_run(self, mesh, mesh_order, omega):
    if self.sol_info['topology']['dim'] == 1:
      base_element = Lobbato1DElement
    else:
      raise ValueError('The dimension is not supported')

    mesh_domains = self.sol_info['topology']['mesh_domain']

    materials = self.sol_info['materials']
    physic_domains = self.sol_info['physic_domain']
    subdomains = {}
    for _, phy_domain in physic_domains.items():
      material = materials[phy_domain[2]]
      domain_elements = mesh_domains[phy_domain[1]]['domain_elements']
      subdomains[material] = domain_elements

    check_material_compability(subdomains)
    elements2node = mesh.mesh_coordinates()

    bases = []
    for mat, elems in subdomains.items():
      if mat.TYPE == 'fluid':
        bases += [
            base_element('Pf', mesh_order[i], elements2node[elem])
            for i, elem in enumerate(elems)
        ]

    fe_space = FESpace(mesh, subdomains, bases)
    Helmholtz_assember = HelmholtzAssembler(fe_space,
                                            subdomains,
                                            dtype=np.complex128)
    Helmholtz_assember.assembly_global_matrix(bases, 'Pf', omega)
    left_hand_matrix = Helmholtz_assember.get_global_matrix()

    right_hand_vec = np.zeros(Helmholtz_assember.nb_global_dofs,
                              dtype=np.complex128)
    BCs_applier = ApplyBoundaryConditions(mesh, fe_space, left_hand_matrix,
                                          right_hand_vec, omega)

    bcs = self.sol_info['boundary_conditions']
    for bc in bcs:
      if bc['type'] == 'fluid_velocity':
        bc['value'] *= np.exp(-1j * omega)
        BCs_applier.apply_nature_bc(bc, var='Pf')
      elif bc['type'] == 'total_displacement':
        BCs_applier.apply_nature_bc(bc, var='Pb')
      elif bc['type'] == 'solid_stress':
        BCs_applier.apply_nature_bc(bc, var='Ux')

    linear_solver = LinearSolver(fe_space=fe_space)
    linear_solver.solve(left_hand_matrix, right_hand_vec)
    sol = linear_solver.u

    return sol

  def post_process(self):
    frf_processer = PostProcessFRF(self.sol_info['frequency'],
                                   r'1D Helmholtz FRF', 'SPL(dB)')
    postprocess_info = self.sol_info['post_processing']
    for post_type, values in postprocess_info.items():
      if post_type == 'frf':
        for frf_id, frf_info in values.items():
          file_name = frf_info['file_name']
          if file_name == 'AUTO':
            file_name = f'frf_{frf_id}.txt'
          frf_processer.save_sol((self.solution_pool, 'simulation'))

      elif post_type == 'map':
        for map_id, map_info in values.items():
          pass
      else:
        raise ValueError('The post processing block is not supported yet')

  def setup(self):
    nodes = self.sol_info['topology']['mesh_nodes']
    connectivity = self.sol_info['topology']['mesh_elements']
    mesh_order = self.sol_info['topology']['mesh_order']
    freqs = self.sol_info['frequency']
    self.solution_pool = defaultdict(dict)
    if self.sol_info['topology']['dim'] == 1:
      mesh = Mesh1D(nodes, connectivity)

    for freq in freqs:
      omega = 2 * np.pi * freq
      lambda_f = self.sol_info['materials'][0].c_f / freq
      h_max = lambda_f / 8
      elem_min_size = mesh.get_min_size()
      while elem_min_size > h_max:
        mesh_order += 1    # p-refinement
        elem_min_size /= mesh_order

      sol = self.fem_run(mesh, mesh_order, omega)
      self.solution_pool[freq] = sol


# wirte the test code for above parser class
# def main():
#     two_fluid_setup = PyAcousiXSetuper('two_fluid.axi')
#     two_fluid_setup.setup()

# main()
