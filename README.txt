---------------------------------------------------------
General Latent Feature Model (GLFM) for matlab and python
---------------------------------------------------------

-----------
Description
-----------

This code implements inference for an Indian Buffet process with heterogeneous
observations. The core code is in C++, but user interfaces in both Matlab and
Python are provided.

You can find some demos of GLFM in action in the demos folder.

------------
Quick start
------------

You can use GLFM from within Matlab or Python.
Below we show an example of how to call GLFM for matrix completion and data
exploration of a given dataset.

Calling from Matlab
-------------------
    bla

Calling from Python
-------------------
    bla

# print help message for required arguments
python -m GLFM.infer --help

------------
Requirements
------------

For Matlab:
    - Matlab 2012b or higher
    - GSL library
        In UBUNTU: sudo apt-get install libgsl0ldbl or sudo apt-get install libgsl0-dev
    - GMP library
        In UBUNTU: sudo apt-get install libgmp3-dev

For Python:
    - Python 2.7 or higher
    - Cython 0.25.2 or higher
    - Libraries: cymsm, cython_gsl

--------------------------
Installation Instructions
--------------------------

To use GLFM for the first time, follow the installation instructions on our
project wiki.
Once installed, please visit the Configuration wiki page to learn how to
configure where data is saved and loaded from on disk.
All documentation can be found on the project wiki.

--------------------------
In order to run GLFM on your data, you need to:

1) Download the latest git repository
2) Install Anaconda, Cython, and a few additional libraries:
    - Anaconda: https://www.continuum.io/downloads
    - Library gsl: conda install gsl
    - Cython: conda install -c anaconda cython=0.25.2
    - Cython_gsl: conda install -c pesoto cython_gsl=1.0.0
    - cymsm library: conda install -c anaconda cymem=1.31.2

3) Compile the C++ code, either for MATLAB or for PYTHON
    - For MATLAB:
        - Add path 'Ccode' and its children to Matlab workspace
        - From matlab command window, execute:
            //mex  -lgsl -lgmp IBPsampler.cpp
            mex  -lgsl -lgmp -lgslcblas IBPsampler.cpp

    - For PYTHON:
        - Go to src/GLFMpython folder
        - run command: python setup.py build_ext --inplace

--------
Citation
--------

Please, cite it as detailed below:
I. Valera, M. F. Pradier and Z. Ghahramani, "Bayesian Nonparametric Latent Feature Model", ArXive, 2016.

-------
Contact
-------

For further information or contact: Isabel Valera at
isabel.valera.martinez@gmail.com or Melanie F. Pradier at melanie.fpradier@gmail.com


