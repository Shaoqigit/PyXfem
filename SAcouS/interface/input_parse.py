from collections import defaultdict
import numpy as np

from SAcouS.acxfem.materials import MaterialFactory

class PyAcousiXSetupParser:
    def __init__(self, file_path):
        self.file_path = file_path
        self.analysis_type = None
        self.frequencies = None
        self.topo_props = {}
        self.dimension = None
        self.mesh_nodes = None
        self.mesh_elements = None
        self.mesh_order = None
        self.domains = {}
        self.materials = {}
        self.physic_domains = []
        self.boundary_conditions = []
        self.solver_type = None
        self.post_processing_requests = []

    def parse_analysis(self):
        pass

    def parse_level1(self):
        blocks = {}
        current_block = None
        with open(self.file_path, 'r') as file:
            for line in file:
                line = line.strip()
                # sklp comments and space lines
                if line.startswith('//') or line == '':
                    continue  # Skip comments

                if line.startswith('# BEGIN'):
                    block_name = line.split('BEGIN ')[1].lower()
                    current_block = []
                    blocks[block_name] = current_block

                    continue

                if line.startswith('# END'):
                    current_block = None
                    continue

                if current_block is not None:
                    current_block.append(line)

            return blocks
     
    def parse_level2plus(self, level:int, blocks:list):
        if level == 2:
            indicator = '##'
        elif level == 3:
            indicator = '###'
        levelp_blocks = {}
        current_block = None
        for line in blocks:
            line = line.strip()
            # sklp comments and space lines
            if line.startswith('//') or line == '':
                continue  # Skip comments

            if line.startswith(indicator+' BEGIN'):
                block_name = line.split('BEGIN ')[1].lower()
                current_block = []
                levelp_blocks[block_name] = current_block
                continue

            if line.startswith(indicator+' END'):
                current_block = None
                continue

            if current_block is not None:
                current_block.append(line)
        return levelp_blocks


    def parse_analysis(self, analysis_blocks: list):
        analysis_block = analysis_blocks[0].split(',')
        self.analysis_type = analysis_block[0]
        self.frequencies = np.linspace(float(analysis_block[1]), float(analysis_block[2]), int(analysis_block[3]))
    
        
    def parse_topo(self, topo_blocks: list):
        for line in topo_blocks:
            if line.startswith('DIMENSION'):
                self.topo_props['dim'] = line.split(',')[1]
                break
        level_2_block = self.parse_level2plus(2, topo_blocks)
        try:
            mesh_block = level_2_block['mesh']
            mesh_block2 = self.parse_level2plus(3, mesh_block)
            # parse mesh nodes
            info_mesh_nodes = mesh_block2['node'][0].split(',')
            if info_mesh_nodes[0] == 'RANGE':
                self.mesh_nodes = np.linspace(float(info_mesh_nodes[1]), float(info_mesh_nodes[2]), int(info_mesh_nodes[3]))
                self.topo_props['mesh_nodes'] = self.mesh_nodes
            elif info_mesh_nodes[0] == 'LIST':
                raise ValueError('The mesh nodes definition has not been implemented yet')
            else:
                raise ValueError('The mesh nodes definition is not correct')
            
            # parse mesh elements
            info_mesh_elements = mesh_block2['element']
            global_order = int(info_mesh_elements[0].split(',')[1])
            mesh_type = info_mesh_elements[1].split(',')[0]
            if mesh_type == 'RANGE':
                start, end, num_nodes = map(float, info_mesh_elements[1].split(','))
                step = (end - start) / (num_nodes - 1)
                self.mesh_elements.append((info_mesh_elements[0].strip(), [start + i * step for i in range(int(num_nodes))]))
            elif mesh_type == 'LIST':
                num_elements = len(info_mesh_elements[2:])
                self.mesh_elements = np.zeros((num_elements, 2), dtype=int)
                self.mesh_order = np.ones((num_elements), dtype=int)
                for i, element in enumerate(info_mesh_elements[2:]):
                    element_info = element.split(',')
                    if element_info[1] == 'NONE':
                        self.mesh_order[i] = global_order
                    else:
                        self.mesh_order[i] = int(element_info[1])
                    self.mesh_elements[i] = np.array([int(node) for node in element_info[2:]])
                self.topo_props['mesh_elements'] = self.mesh_elements
                self.topo_props['mesh_order'] = self.mesh_order
        except KeyError:
            print('No mesh block is defined')
        
        try:
            domain_block = level_2_block['domain']
            self.topo_props[f'mesh_domain'] = list()
            for domain in domain_block:
                domain_info = domain.split(',')
                domain_id = domain_info[0]
                domain_name = domain_info[1]
                domain_elements = np.array([int(element.strip()) for element in domain_info[2:]])
                self.topo_props['mesh_domain'].append({'domain_id': domain_id, 'domain_,name': domain_name, 'domain_elements': domain_elements})
        except KeyError:
            print('No domain block is defined')


    def parse_materials(self, material_blocks: list):
        build_material = MaterialFactory()
        for material in material_blocks:
            material_info = material.split(',')
            material_id = int(material_info[0])
            material_type = material_info[1]
            material_name = material_info[2]
            mat_properties = material_info[3:]
            if material_type == 'AIR':
                properties_values = []
            else:
                properties_values = [float(prop.strip()) for prop in mat_properties]
            self.materials[material_id] = build_material.create_material(material_type, material_name, *properties_values)


    def parse_physic_domains(self, physic_domain_blocks: list):
        for physic_domain in physic_domain_blocks:
            physic_type, domain_id, material_id = physic_domain.split(',')
            self.physic_domains.append((physic_type.strip(), int(material_id.strip()), int(domain_id.strip())))
        

    def parse_boundary_conditions(self, bc_blocks: list):
        for bc in bc_blocks:
            bc_info = bc.split(',')
            bc_type = bc_info[0]
            bc_domain = bc_info[1]
            bc_value = bc_info[2]
            self.boundary_conditions.append((bc_type.strip(), int(bc_domain.strip()), float(bc_value.strip())))

    def parse_solver(self, line):
        _, solver_type = line.split(',')
        self.solver_type = solver_type.strip()

    def parse_post_processing_requests(self, line):
        _, request_type = line.split(',')
        self.post_processing_requests.append(request_type.strip())

    
# wirte the test code for above parser class
def test_parser():
    parser = PyAcousiXSetupParser('two_fluid.axi')
    blocks = parser.parse_level1()
    parser.parse_analysis(blocks['analysis'])
    parser.parse_topo(blocks['topology'])
    parser.parse_materials(blocks['material'])

test_parser()

