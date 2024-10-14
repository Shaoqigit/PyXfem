from abc import abstractmethod, ABCMeta
import numpy as np

from SAcouS.Materials import MaterialFactory

class BaseParser(metaclass=ABCMeta):
    def __init__(self, file_path):
        self.file_path = file_path

    @abstractmethod
    def parse(self):
        pass


class AcoustiXPaser(BaseParser):
    def __init__(self, file_path):
        super().__init__(file_path)
        if not file_path.endswith('.axi'):
            raise ValueError('File must be an .axi file')
        self.analysis_type = None
        self.frequencies = None
        self.topology_props = {}
        self.materials = {}
        self.physic_domains = {}
        self.boundary_conditions = {}
        self.solver = {}
        self.post_processing = {}

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

    def parse_topology(self, topo_blocks: list):
        for line in topo_blocks:
            if line.startswith('DIMENSION'):
                self.topology_props['dim'] = int(line.split(',')[1])
                break
        level_2_block = self.parse_level2plus(2, topo_blocks)
        try:
            mesh_block = level_2_block['mesh']
            mesh_block2 = self.parse_level2plus(3, mesh_block)
            # parse mesh nodes
            info_mesh_nodes = mesh_block2['node'][0].split(',')
            if info_mesh_nodes[0] == 'RANGE':
                self.mesh_nodes = np.linspace(float(info_mesh_nodes[1]), float(info_mesh_nodes[2]), int(info_mesh_nodes[3]))
                self.topology_props['mesh_nodes'] = self.mesh_nodes
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
                self.topology_props['mesh_elements'] = self.mesh_elements
                self.topology_props['mesh_order'] = self.mesh_order
        except KeyError:
            print('No mesh block is defined')
        
        try:
            domain_block = level_2_block['domain']
            self.topology_props['mesh_domain'] = {}
            for domain in domain_block:
                domain_info = domain.split(',')
                domain_id = int(domain_info[0])
                domain_dim = domain_info[1]
                domain_name = domain_info[2]
                domain_elements = np.array([int(element.strip()) for element in domain_info[3:]])

                self.topology_props['mesh_domain'][domain_id] = {'domain_dim': domain_dim, 'domain_,name': domain_name, 'domain_elements': domain_elements}
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
        for i, physic_domain in enumerate(physic_domain_blocks):
            physic_type, domain_id, material_id = physic_domain.split(',')
            self.physic_domains[i] = [physic_type.strip(), int(material_id.strip()), int(domain_id.strip())]

    def parse_boundary_conditions(self, bc_blocks: list):
        for bc in bc_blocks:
            bc_info = bc.split(',')
            bc_id = int(bc_info[0])
            bc_type = bc_info[1].lower()
            bc_domain = int(bc_info[2].strip())
            bc_value = float(bc_info[3].strip())
            self.boundary_conditions[bc_id]={'type':bc_type.strip(), 'position':bc_domain, 'value': bc_value}

    def parse_solver(self, solver_blocks: list):
        for solver_block in solver_blocks:
            solver_info = solver_block.split(',')
            solver_id = int(solver_info[0])
            solver_type = solver_info[1]
        self.solver[solver_id] = solver_type

    def parse_post_processing(self, post_pro_blocks: list):
        frfs_out = {}
        maps_out = {}
        level_2_block = self.parse_level2plus(2, post_pro_blocks)
        for block_name, block in level_2_block.items():
            if block_name == 'frf':
                for frf in block:
                    frf_info = frf.split(',')
                    frf_id = int(frf_info[0])
                    frfs_out[frf_id]={'file_name': frf_info[1].strip(), 'domain': int(frf_info[2].strip())}
            elif block_name == 'map':
                for map_out in block:
                    map_out_info = map_out.split(',')
                    map_id = int(map_out_info[0])
                    maps_out[map_id]={'file_name': map_out_info[1].strip(), 'freqs': [float(freq) for freq in map_out_info[2:]]}
            else:
                raise ValueError('The post processing block is not supported yet')
        self.post_processing = {'frfs': frfs_out, 'maps': maps_out}


    def parse(self):
        blocks = self.parse_level1()
        for bock_name, block in blocks.items():
            if bock_name == 'analysis':
                self.parse_analysis(block)
            elif bock_name == 'topology':
                self.parse_topology(block)
            elif bock_name == 'material':
                self.parse_materials(block)
            elif bock_name == 'physic_domain':
                self.parse_physic_domains(block)
            elif bock_name == 'boundary_condition':
                self.parse_boundary_conditions(block)
            elif bock_name == 'solver':
                self.parse_solver(block)
            elif bock_name == 'post_pro':
                self.parse_post_processing(block)
            else:
                raise ValueError('The block name is not supported')
            
        return {'analysis': self.analysis_type, 'frequencies': self.frequencies, 'topology': self.topology_props, 'materials': self.materials, 'physic_domain': self.physic_domains, 'BCs': self.boundary_conditions, 'solver': self.solver, 'post_processing': self.post_processing}
    

class ParserFactory:

    def __init__(self):
        pass

    @staticmethod
    def create_parser(file_path):
        # import pdb;pdb.set_trace()
        if file_path.endswith('.axi'):
            return AcoustiXPaser(file_path)
        else:
            raise ValueError('other types of input files are not supported yet')