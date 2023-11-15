import sys
sys.path.append('../src')
import os
import subprocess



test_path = "tests/"
test_cases = ['test_material_pem.py', 
              'test_absorption_comp.py',
              'test1_two_layer.py', 
              'test_two_fluid_new.py', 
              'test2_impedance_bc.py', 
              'test_biot_equation.py',
              'test_biot_equation_new.py',
              'test_modal_reduction_FRF.py',
              'test_fadm.py',
              'test_tmm_biot.py']
for test_case in test_cases:
    cmd = ['python3', test_path+test_case]
    # print(cmd)
    result = subprocess.run(cmd, stdout=subprocess.PIPE)
    if 'Test passed!' in result.stdout.decode('utf-8'):
        # import pdb; pdb.set_trace()
        print(f"Test case:  {test_case:<40}", f"\033[1;32m {'SUCCESS':>30}\033[0m")
    else:
        print(f"Test case:  {test_case:<40}", f"\033[1;31m {'FAILED' :>30}\033[0m")
