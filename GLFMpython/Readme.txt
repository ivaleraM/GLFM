# ----------------------------------------------------
Python Wrapper for General Latent Feature Model (GLFM)
# ----------------------------------------------------

------------
Requirements
------------

You should install Cython and the C++ GSL library in order to be able to call these functions.

-----------
Compilation
-----------

Let us first generate the python wrapper, gsl_run.cpp
For that, go to ../Ccode/wrapper_python/
Run command:
    python setup.py build_ext --inplace

A new library object gsl_run.so will be generated. Now you can import such
library in your python scripts in order to run the GLFM model.

See demo scripts in this folder:

    main_toy.py
    simMiss.py
    matrix_completion.py

