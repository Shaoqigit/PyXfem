# **S**imple **A**coustic **S**imulator(SAcouS)
[comment]: <> ( - AcoustiX based Sound simulation library)
[comment]: <> ( - put following badges centre)
<center>
  
[![PyPi Version](https://img.shields.io/pypi/v/pyacoustix.svg?style=flat-square)](https://pypi.org/project/PyAcoustiX/)

</center>

                                 based on Finite Element Method
                                 contact: Shaoqiwu@outlook.com

**SAcouS** is pure python software based on an Finite Element libraray 'PyacoustiX' (AcoustiX in cpp version found in my git repo), mainly designed for Acoustic simulation but can be used (extended) to other physics. This labrary aims at providing the simplest, clearest and easily conprehensive tool to use FEM for purpose of research, education or industry prototype.

## Installation and Usage
```bash
pip install pyacoustix
```
- You can use either as a standalone library or as a part of your project. The following is a simple example to use the library:
```python
from SAcouS import *
```
- either used with executable with standard pyacoustix input file format .axi. To do so, you need first install the package as
```bash
pip install .
```
Then you can run the following command to execute the simulation:
```bash
sacous -i input_file.axi
```
- This software is also a wrap of python interpreter, so you can use it interactively with the comsole by entering:
```bash
sacous
```
you get into PyacoustiX console, where you can use the library interactively. For example:
```bash
Welcome to PyAcoustiX console!
PyAcoustiX 0.9.9
base on python 3.10.13 (main, Sep 11 2023, 13:44:35) [GCC 11.2.0] on linux
Type "help", "copyright" or "license" for more information.
Author: Shaoqiwu@outlook.com
(PyAcoustiXInteractiveConsole)
>>> import numpy as np
```

## Main special features/architecture of the library compared to other existing lib:
> * Apart from classical Linear/Quadratic Larange polynomial, **High order (lobatto)** shape functions are supported.
> * The element type (interpolation order of shape functions) can be **varied** on each element (mesh).
> * **Extended elements** and double node techinique are used to deal with discontitnuity.
> * Various common Acoustic materials are supported, including porous rigid, Limp and **Biot-up** models.
> * FEM-**TMM** coupling is supported to model multi-layered configuration in general geometry problems.


### Up to this moment:
* High order lobatto shape function supported (up to $p=4$)
* JCA and Limp model to account for porous acoustic materials
* **Biot equation** modeling for poro-elastic materials
* Impedence boundary condition supported
* Weakly enforced essential boundary condition (penalty method and **Nitche's method**)
* Biot UP and Fluid coupling model (implemented but not fully validated)
* Modal domain reduction method

### To do list
* update to python 3.11
* MoR for damped systems (RB+EIM algorithme for porous materials)
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
Test case:  main_test4_biot.py                                           SUCCESS
Test case:  main_test5_biot_up.py                                        SUCCESS
Test case:  main_absorption_comp.py                                      SUCCESS
```
