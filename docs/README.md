[MATLAB functions](doc_matlab.html) | [PYTHON functions](doc_matlab.html) | [FAQ](FAQ_errors.html)

Description
-----------

This code implements inference for an Indian Buffet process with heterogeneous
observations. The core code is in C++. User interfaces for both MATLAB and PYTHON are provided.

You can find some demos of GLFM in action in the demos folder.

Quick start
------------

You can use GLFM from within Matlab or Python.
Below we show an example of how to call GLFM for matrix completion and data
exploration of a given dataset.

**Calling from Matlab**

    hidden = IBPsampler_run(data);

where data is a structure containing:
    X: ´N x D´ observation matrix of N samples and D dimensions
    C: ´1 x D´ string array indicating type of data for each dimension

--- Alternative calls ---

    hidden = IBPsampler_run(data, hidden);
OR
    hidden = IBPsampler_run(data, hidden, params);

where hidden is a structure of latent variables:
    Z: ´N xK´ binary matrix of feature assignments (initialization for the IBP)
and params is a structure containing all simulation parameters and model
    hyperparameters (see [description of structures](doc_struct.html) for further details).

Check further documentation for available [MATLAB functions](doc_matlab.html).

**Calling from Python**

    import GLFM
    (Z_out,B_out,Theta_out) = GLFM.infer(X,C,Z)

Several optional parameters can be provided in the function call, check further documentation for available [PYTHON functions](doc_python.html).

Requirements
-------------

* For Matlab:
    * Matlab 2012b or higher
    * MEX environment (to compile C++ from Matlab terminal)
    * C++ compiler (typically already installed)
    * Libraries: gsl, gmp

* For Python:
    * Python 2.7 or higher
    * Cython 0.25.2 or higher
    * C++ compiler (typically already installed)
    * Libraries: cymsm, cython_gsl

Installation Instructions
--------------------------

In order to run GLFM on your data, you need to:

1. Download the latest git repository
2. Install the necessary libraries:

    * **For Python**
        * Anaconda: https://www.continuum.io/downloads
        * Cython: http://cython.org
        * Libraries gsl, cython_gsl and cymsm

    Once Anaconda is available, you can directly install all the rest as:

        conda install -c anaconda cython=0.25.2
        conda install gsl
        conda install -c pesoto cython_gsl=1.0.0
        conda install -c anaconda cymem=1.31.2

    * **For Matlab**
        * GSL library
        * GMP library

    In Ubuntu, these libraries can be easily installed executing from a terminal:

        sudo apt-get install libgsl0ldbl
        sudo apt-get install libgmp3-dev
            OR
        sudo apt-get install libgsl0-dev
        sudo apt-get install libgmp3-dev

3. Compile the C++ code, either for MATLAB or for PYTHON

For MATLAB, you should add the path 'Ccode' and its children to Matlab workspace. Then, from matlab command window, execute:

    mex  -lgsl -lgmp -lgslcblas IBPsampler.cpp

For PYTHON, go to src/GLFMpython folder and execute the following command from terminal:

    python setup.py build_ext --inplace

Check documented list of [problems at installation time](FAQ_errors.md) and how to solve them.

Citation and contact
--------------------

Please, cite it as detailed below:
> I. Valera, M. F. Pradier and Z. Ghahramani, "Bayesian Nonparametric Latent Feature Model", ArXive, 2016.

For further information or contact:
* Isabel Valera: isabel.valera.martinez (at) gmail.com
* Melanie F. Pradier: melanie.fpradier (at) gmail.com


