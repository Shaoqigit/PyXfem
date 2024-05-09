from collections import defaultdict
import numpy as np
class PyAcousiXSetupParser:
    def __init__(self, file_path):
        self.file_path = file_path
        self.analysis_type = None
        self.frequency_range = None
        self.dimension = None
        self.mesh_nodes = []
        self.mesh_elements = []
        self.domains = []
        self.materials = []
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
                # import pdb; pdb.set_trace()
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
        blocks = {}
        current_block = None
        for line in blocks:
            # import pdb; pdb.set_trace()
            line = line.strip()
            # sklp comments and space lines
            if line.startswith('//') or line == '':
                continue  # Skip comments

            if line.startswith(indicator+' BEGIN'):
                block_name = line.split('BEGIN ')[1].lower()
                current_block = []
                blocks[block_name] = current_block
                continue

            if line.startswith(indicator+' END'):
                current_block = None
                continue

            if current_block is not None:
                current_block.append(line)
        return blocks


    def parse_analysis(self, analysis_blocks: list):
        analysis_block = analysis_blocks[0].split(',')
        analysis_type = analysis_block[0]
        frequency_range = np.linspace(analysis_block[1], analysis_block[2], analysis_block[3])
        return analysis_type, frequency_range
    
        
    def parse_topo(self, topo_blocks: list):
        topo_props = {}
        blocks = None
        current_block = None
        for line in topo_blocks:
            if line.startswith('DIMENSION'):
                topo_props['dim'] = line.split(',')[1]
                break
        self.parse_level2plus(2, topo_blocks)

    def parse_dimension(self, line):
        _, dimension = line.split(',')
        self.dimension = int(dimension.strip())

    def parse_mesh_nodes(self, line):
        _, node_range = line.split(',')
        start, end, num_nodes = map(float, node_range.split(','))
        step = (end - start) / (num_nodes - 1)
        self.mesh_nodes = [start + i * step for i in range(int(num_nodes))]

    def parse_mesh_elements(self, line):
        _, element_type, *elements = line.split(',')
        elements = [int(element.strip()) for element in elements]
        self.mesh_elements.append((element_type.strip(), elements))

    def parse_domains(self, line):
        _, domain_name, *element_ids = line.split(',')
        element_ids = [int(element_id.strip()) for element_id in element_ids]
        self.domains.append((domain_name.strip(), element_ids))

    def parse_materials(self, line):
        pass
        # _, material_type, material_name = line.split(',')
        # properties = {}
        # for line in lines:
        #     line = line.strip()
        #     if line.startswith('# END MATERIAL'):
        #         break
        #     elif line.startswith('//') or line.startswith('##') or line.startswith('###'):
        #         continue
        #     else:
        #         property_name, property_value = line.split(',')
        #         properties[property_name.strip()] = property_value.strip()
        # self.materials.append((material_type.strip(), material_name.strip(), properties))

    def parse_physic_domains(self, line):
        _, physic_type, material_id, domain_id = line.split(',')
        self.physic_domains.append((physic_type.strip(), int(material_id.strip()), int(domain_id.strip())))

    def parse_boundary_conditions(self, line):
        _, boundary_type, value, domain_id = line.split(',')
        self.boundary_conditions.append((boundary_type.strip(), float(value.strip()), int(domain_id.strip())))

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
    parser.parse_topo(blocks['analysis'])
    import pdb; pdb.set_trace()

test_parser()

