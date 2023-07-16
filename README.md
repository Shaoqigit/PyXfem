### SAcSim: Simple Acoustic Simulator
                   based on Finite Element Method
                    contact: Shaoqiwu@outlook.com

SAcSim is Finite Element (FEM) modular Library, mainly designed for Acoustic simulation but can be used (extended) to other physics. This labrary aims at providing the simplest, clearest and easily conprehensive tool to use FEM for purpose of research, education or industry prototype.

Main special features/architecture of the library compared to other existing lib:
* Apart from classical Linear/Quadratic Larange polynomial, High order (lobatto) shape functions are supported
* The element type (interpolation order of shape functions) can be variable on each element (mesh)
* Extended elements and double node techinique are used to deal with discontitnuity
* Various common Acoustic materials are supported, including porous rigid, Limp and Biot-up models.

Up to this moment:

* High order lobatto shape function supported (up to $p=4$)
* JCA and Limp model to account for porous acoustic materials
* Impedence boundary condition supported

Roadmap:

* Biot UP solid and Fluid coupling model
* Perfect Matched Layer (PML) for free field boundary condition
* Infinite Element for free field boundary condition
