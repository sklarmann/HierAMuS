# HierAMuS - A hierarchical higher order and multiscale simulation finite element program

I started this project around 2016 as a Ph.D. student as a private project to learn C++. I dropped the program for some years. In the meantime, I started reusing it, realizing that it offers great flexibility in designing different finite element structures.

As the code matured, I decided to make it publicly available. It still has a stage of an early alpha version. As it was a learning project, many things are not implemented in a suitable way.
Therefore, expect a lot of refactoring of the code in the future.


# Installation
The program is managed by CMake. It requires an installtion of Python with the package pybind11. 
During the configuration process, CMake will download the necessary libraries, compiles them and installs them in the source folder under stage. 
As one dependent library is Paraview, the first configuration process in CMake might take around 1 hour.

After creating the projekt with CMake, the build target PyPackageInstall installs the package as an development version in the current python environment. 
When using linux, this build project ceates an virutal environment in the source folder.

During the python installation the package matplotlib, gmsh and numpy will be installed in the python environment if they are not present.

# Required third-party libraries
- Eigen >= 3.4
- Spectra
- Paraview
- spdlog
- pybind11


# Current publications based on this code

- Coupling 2D continuum and beam elements: a mixed formulation for avoiding spurious stresses DOI: 10.1007/s00466-022-02221-7
- Possibilities and drawbacks using arbitrary precision numbers forstructural analysis DOI: 10.1002/pamm.202000079 (Support of arbitrary precision is currently removed)

# Further comments
In the mean time, more people started to contribute. 
I want to honor their effort and they are listet in the CONTRIBUTORS.md