General Latent Feature Model (GLFM)

Compile
Before using the code, please add path 'C\ code' and compile the sampler function in Matlab with the command:
mex  -lgsl -lgmp IBPsampler.cpp

Table completion function
The function matrix_completion.m completes the table with the MAP solution for the missing data. Please refer to the function for details on the inputs and outputs of the function.

Reproduce experiments in the paper 
The script simMiss.m reproduce the experiment on the QSAR biodegradation database.

Description
Please refer to Description.m

Citation
Please, cite it as detailed below:
I. Valera and Z. Ghahramani, "General Table Completion using a Bayesian Nonparametric Model", Neural Information Processing Systems Conference 2014 (NIPS 2014). Montreal (Canada), 2014.

Contact
For further information or contact Isabel Valera at ivalera@mpi-sws.org.