# AcoustiX Simple Acoustic Simulator
                   based on Finite Element Method
                    contact: Shaoqiwu@outlook.com

**AcoustiX** is Finite Element (FEM) modular Library, mainly designed for Acoustic simulation but can be used (extended) to other physics. This labrary aims at providing the simplest, clearest and easily conprehensive tool to use FEM for purpose of research, education or industry prototype.

## Main special features/architecture of the library compared to other existing lib:
* Apart from classical Linear/Quadratic Larange polynomial, **High order (lobatto)** shape functions are supported
* The element type (interpolation order of shape functions) can be **varied** on each element (mesh)
* **Extended elements** and double node techinique are used to deal with discontitnuity
* Various common Acoustic materials are supported, including porous rigid, Limp and **Biot-up** models.
* FEM-**TMM** coupling is supported to model multi-layered configuration in general geometry problems
* The proposed Finite Admittance Method **(FAM)** (experimental) is supported

### Up to this moment:
* High order lobatto shape function supported (up to $p=4$)
* JCA and Limp model to account for porous acoustic materials
* **Biot equation** modeling for poro-elastic materials
* Impedence boundary condition supported
* Weakly enforced essential boundary condition (penalty method and **Nitche's method**)
* Biot UP and Fluid coupling model (implemented but not fully validated)

### To do list
* Perfect Matched Layer (PML) for free field boundary condition
* Infinite Element for free field boundary condition
* Finite Admittance Method (FAM) for high precision multi-layered configuration
* TMM embeded for multi-layered configuration
* higher order shape function and integration ($p=10$)

### Roadmap:
* Message and log system

## For development and contribution
Every time one new feature is developed, a new test case has to be added in the /tests folder. The run_test must to be run and passed with all green light before created Pull Request as:
```bash
Test case:  main_test1_two_layer.py                                      SUCCESS
Test case:  main_test2_impedance_bc.py                                   SUCCESS
Test case:  main_test3_lobatto.py                                        SUCCESS
```
