import numpy as np
import meshio
from src.basis import Lobbato1DElement
from src.mesh import Mesh1D, DofHandler1D
from src.assembly import Assembler


nodes = np.linspace(-1, 1, 5)

elements = np.array([[0, 1],
                     [1, 2],
                     [2, 3],
                     [3, 4]])

# read the mesh data structure
mesh = Mesh1D(nodes, elements)
mesh.read_mesh()
elements_set = mesh.mesh  # elements with nodes coodinates
print(mesh.mesh)

bases = []  # basis applied on each element, could be different order and type
order = 1
# assemble the basis on the elements
for key, elem in elements_set.items():
    print(elem)
    basis = Lobbato1DElement(order, elem)
    print("local_dofs_index", basis.local_dofs_index())
    print("basis internam dof", basis.num_internal_dofs())
    bases.append(basis)
    order += 1

# handler the dofs
dof_handler = DofHandler1D(mesh, bases)
print(dof_handler.global_dof_index())
print(dof_handler.get_num_dofs())

# initialize the assembler
assemble = Assembler(dof_handler, bases)

assemble.assemble_K()

print(assemble.get_matrix_in_array())
# print(assemble.M)

