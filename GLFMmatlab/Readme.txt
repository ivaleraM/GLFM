General Latent Feature Model (GLFM)

Requirements
In order to compile, it is necessary the C++ library GSL.

Compile
Before using the code, please add path 'Ccode' and compile the sampler function in Matlab with the command:
mex  -lgsl -lgmp IBPsampler.cpp

Requirements:
You need to install first the libraries gsl and gmp in your computer:

UBUNTU:
    sudo apt-get install libgsl0ldbl or sudo apt-get install libgsl0-dev
    sudo apt-get install libgmp3-dev

----- Scripts --------
IBPsampler.cpp
Implements in C, the inference for the IBP with heterogeneous observations. For the description of the inputs and outputs, please refer to Description.m

matrix_completion.m
The function matrix_completion.m completes the table with the MAP solution for the missing data. Please refer to the function for details on the inputs and outputs of the function.

simMiss.m
Reproduce experiments of the paper on Matrix completion on the QSAR biodegradation database.



