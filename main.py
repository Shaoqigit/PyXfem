import numpy as np
import meshio
import matplotlib.pyplot as plt
from matplotlib.pyplot import spy

from src.basis import Lobbato1DElement
from src.mesh import Mesh1D
from src.dofhandler import DofHandler1D
from src.assembly import Assembler
from src.materials import Air, Fluid, EquivalentFluid, check_material_compability


nodes = np.linspace(-1, 1, 5)

elements = np.array([[0, 1],
                     [1, 2],
                     [2, 3],
                     [3, 4]])



# read the mesh data structure
mesh = Mesh1D(nodes, elements)
mesh.refine_mesh(2)

elements_set = mesh.read_mesh()  # elements with nodes coodinates
print(elements_set)

bases = []  # basis applied on each element, could be different order and type
order = 2
# assemble the basis on the elements
for key, elem in elements_set.items():
    print(elem)
    basis = Lobbato1DElement(order, elem)
    print("local_dofs_index", basis.local_dofs_index())
    print("basis internam dof", basis.num_internal_dofs())
    bases.append(basis)

# handler the dofs
dof_handler = DofHandler1D(mesh, bases)
print(dof_handler.global_dof_index())
print(dof_handler.get_num_dofs())
print(dof_handler.num_internal_dofs())
print(dof_handler.num_external_dofs())


# define the materials
air = Air('classical air')
water=Fluid('water', 997, 1481)
phi          = 0.98  # porosity
sigma        = 3.75e3 # resistivity
alpha        = 1.17  # Tortuosity
Lambda_prime = 742e-6  # Viscous characteristic length
Lambda      = 110e-6  # 
xfm = EquivalentFluid('xfm', phi, sigma, alpha, Lambda_prime, Lambda)

# define the subdomains: domain name (material) and the elements in the domain
subdomains = {air: [0,1,2,3,4,5,6,7,8,9], water: [10,11,12,13,14,15]}
check_material_compability(subdomains)


# initialize the assembler
assembler = Assembler(dof_handler, bases, np.float64)


freq = 1000
omega = 2 * np.pi * freq
K_g= assembler.assemble_material_K(subdomains, omega)

M_g= assembler.assemble_material_M(subdomains, omega)
# assemble.assemble_K()
# print(K_g)
spy(assembler.get_matrix_in_array(K_g))# print(assemble.M)
plt.show()

# construct linear system
left_hand_matrix = K_g-M_g
#  natural boundary condition   
nature_bcs = {'type': 'velocity', 'value': 1, 'position': 0}
right_hand_vector = assembler.assemble_nature_bc(nature_bcs)

