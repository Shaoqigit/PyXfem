import sys

sys.path.append('../src')
import os
import subprocess
import concurrent.futures

test_path = "tests/"
test_cases = [
    'test_2D_matrices.py', 'test_2D_air.py', 'test_material_pem.py',
    'test_absorption_comp.py', 'test1_two_layer.py', 'test_two_fluid_new.py',
    'test2_impedance_bc.py', 'test_biot_equation.py',
    'test_biot_equation_new.py', 'test_modal_reduction_FRF.py', 'test_fadm.py',
    'test_tmm_biot.py'
]
# remove result files

import glob
for file in glob.glob(test_path + '*.pos'):
  os.remove(file)
# search in the directory to see if failed_test_cases.txt exists
# if it exists, reead it and run the failed test cases
# if not, run all the test cases
# if the test case is passed, print the test case name in green
# if the test case is failed, print the test case name in red
faield_test_cases = []
if os.path.exists(test_path + 'failed_test_cases.txt'):
  with open(test_path + 'failed_test_cases.txt', 'r') as f:
    for line in f:
      faield_test_cases.append(line.strip())

if faield_test_cases:
  test_cases = faield_test_cases

faield_test_cases = []


# run the test cases parallelly
def run_test_case(test_case):
  cmd = ['python3', test_path + test_case]
  result = subprocess.run(cmd, stdout=subprocess.PIPE)
  if 'Test passed!' in result.stdout.decode('utf-8'):
    print(f"Test case:  {test_case:<40}", f"\033[1;32m {'SUCCESS':>30}\033[0m")
  else:
    print(f"Test case:  {test_case:<40}", f"\033[1;31m {'FAILED' :>30}\033[0m")
    return test_case


# ProcessPoolExecutor
# ThreadPoolExecutor
import time

start_time = time.time()
with concurrent.futures.ProcessPoolExecutor() as executor:
  results = executor.map(run_test_case, test_cases)

faield_test_cases = [result for result in results if result is not None]

if not faield_test_cases:
  print(f"\033[1;32m {'All test cases passed!':>30}\033[0m")
  try:
    os.remove(test_path + 'failed_test_cases.txt')
  except FileNotFoundError:
    pass
else:
  with open(test_path + 'failed_test_cases.txt', 'w') as f:
    for test_case in faield_test_cases:
      f.write(test_case + '\n')

print(f"Time taken: {time.time() - start_time:.2f} seconds")
