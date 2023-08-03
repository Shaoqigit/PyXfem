import sys
sys.path.append('../src')
import os
import subprocess



test_path = "tests/"
test_cases = ['test_material_pem.py', 
              'main_test1_two_layer.py', 
              'main_test2_impedance_bc.py', 
              'main_test_biot_equation.py']
for test_case in test_cases:
    cmd = ['python3', test_path+test_case]
    print(cmd)
    result = subprocess.run(cmd, stdout=subprocess.PIPE)
    if 'Test passed!' in result.stdout.decode('utf-8'):
        # import pdb; pdb.set_trace()
        print(f"Test case:  {test_case:<40}", f"\033[1;32m {'SUCCESS':>30}\033[0m")
    else:
        print(f"Test case:  {test_case:<40}", f"\033[1;31m {'FAILED' :>30}\033[0m")
