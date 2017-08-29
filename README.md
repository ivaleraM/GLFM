GLFM: General Latent Feature Modeling toolbox for matlab and python
---------------------------------------------------------

This code implements a package for General Laten Feature Model (GLFM) suitable for heterogeneous
observations. The core code is in C++. User interfaces in both Matlab and
Python are provided. Moreover, several demos are provided to illustrate different applications, including missing data estimation and data exploratory analysis, of the GLFM.

To cite this work, please use

     I. Valera, M. F. Pradier and Z. Ghahramani, 
     "General Latent Feature Model for Heterogeneous Datasets", 
     Available on [ArXiv](https://arxiv.org/abs/1706.03779), 2017.

-----------------
GLFM Description
-----------------

GLFM is a general Bayesian nonparametric latent feature model suitable for heterogeneous datasets, where the attributes describing each object can be either discrete, continuous or mixed variables. Specifically, it accounts for the following types of data:

• Continuous variables:

* Real-valued (encoded as 'g'): the attribute takes values in the real line. 	
* Positive real-valued ('p'): the attribute takes values in the real line.

• Discrete variables:
* Categorical data ('c'): the attribute takes a value in a finite unordered set, e.g., {‘blue’,‘red’, ‘black’}.
* Ordinal data ('o'): the attribute takes values in a finite ordered set, e.g., {‘never’, ‘often’, ‘always’}.
* Count data ('n'): the attribute takes values in the set {0,...,∞}.

The GLFM builds on the Indian Buffet Process (Griffiths and Ghahramani, 2011), and therefore, it assumes that each observation x_n^d can be explained by  a potentially infinite-length binary vector **z**_n whose elements indicate whether a latent feature is active or not for the n-th object; and a (real-valued) weighting vector **B**^d, whose elements weight the influence of each latent feature in the d-th attribute. 
Since the product of the latent feature vector and the weighting vector leads to a real-valued variable, it is necessary to map this variable to the desirable output (continuous or discrete) space, for example, the positive real line. Thus, the GLFM assumes the existence of intermediate Gaussian variables y_n^d, with mean **z**_n**B**^d and called pseudo-observation, and a transformation function f_d() that maps this variable into the actual observation x_n^d.

As an example, an ordinal attribute taking values in the ordered set {*low*, *medium*, *high*} can be represented using the GLFM as:


<p align="center">
  <img src="https://raw.githubusercontent.com/ivaleraM/GLFM/master/figures/Model_example.png" width="600"/>
  <!-- <img src="https://github.com/ivaleraM/GLFM/blob/master/figures/Model_example.png" width="600"/> -->
</p>

For more details on the GLMF, please refer to the paper. 

------------
GLFM Toolbox
------------

You can use GLFM from within Matlab or Python. Below we show an example of how to run the GLFM inference. For mode details on the functions and data structures, please refer to the [GLFM documentation](https://ivaleram.github.io/GLFM/).


Calling from Python
-------------------
    import GLFM
    (hidden) = GLFM.infer(data)

where **data** is a structure containing:

    X: NxD observation matrix of N samples and D dimensions

    C: 1xD string array indicating type of data for each dimension 


--- Alternative calls ---
    
    import GLFM
    hidden = GLFM.infer(data, hidden);
OR

    import GLFM
    hidden = GLFM.infer(data, hidden, params);

where **hidden** is a structure of latent variables:

    Z: NxK binary matrix of feature assignments (initialization for the IBP)

and **params** is a structure containing all simulation parameters and model hyperparameters (see documentation for further details).

Calling from Matlab
-------------------
    hidden = GLFM_infer(data);

where **data** is a structure containing:

    X: NxD observation matrix of N samples and D dimensions

    C: 1xD string array indicating type of data for each dimension 



--- Alternative calls ---

    hidden = GLFM_infer(data, hidden);
OR

    hidden = GLFM_infer(data, hidden, params);

where **hidden** is a structure of latent variables:

    Z: NxK binary matrix of feature assignments (initialization for the IBP)

and **params** is a structure containing all simulation parameters and model
    hyperparameters (see documentation for further details).

Calling from R
---------------
    hidden <- GLFM_infer(data)

where **data** is a structure containing:

    X: NxD observation matrix of N samples and D dimensions

    C: 1xD string array indicating type of data for each dimension 



--- Alternative calls ---

    hidden <- GLFM_infer(data,hidden)
OR

    hidden = GLFM_infer(data, list(hidden, params));

where **hidden** is a list of latent variables:

    Z: NxK binary matrix of feature assignments (initialization for the IBP)

and **params** is a list containing all simulation parameters and model
    hyperparameters (see documentation for further details).


Requirements
------------

For Python:

    - Python 2.7 or higher
    - Cython 0.25.2 or higher
    - Libraries: cymsm, cython_gsl

For Matlab:

    - Matlab 2012b or higher
    - GNU GSLlibrary
        In UBUNTU: sudo apt-get install libgsl0ldbl or sudo apt-get install libgsl0-dev
    - GMP library
        In UBUNTU: sudo apt-get install libgmp3-dev

For R:

    - R or Rstudio
    - GNU GSL library (e.g. libgsl0-dev on Debian or Ubuntu)
    - Rcpp for seamless R and C++ integration


Installation Instructions
--------------------------

In order to run GLFM on your data, you need to:

1) Download the latest git repository

2) Install necessary libraries:

    - For PYTHON:
    	- Anaconda: https://www.continuum.io/downloads
    	- Library gsl: conda install gsl
    	- Cython: conda install -c anaconda cython=0.25.2
    	- Cython_gsl: conda install -c pesoto cython_gsl=1.0.0
    	- cymsm library: conda install -c anaconda cymem=1.31.2

3) Compile the C++ code as

    - For PYTHON In a terminal):
        - Go to GLFM/src/Ccode/wrapper_python/ folder
        - Run command: python setup.py build_ext --inplace

     - For MATLAB:
        - Add path 'GLFM/src/Ccode' and its children to Matlab workspace
        - From matlab command window, execute: mex  -lgsl -lgmp -lgslcblas IBPsampler.cpp

     - For R (in R or Rstudio):
        - Set directory as: >> setwd("GLFM/src/Ccode")
        - Compile using the following command liens:
            >> install.packages("GLFM/src/Ccode/wrapper_R", repos = NULL, type="source") 
            >> require("Rcpp")
            >> compileAttributes("GLFM/src/Ccode/wrapper_R",verbose=TRUE)
        - Set the working directory as: >> setwd("~/your directory/GLFM/src/GLFMR")


------------
GLFM Demos
------------
The folder `demos' contain scripts, as well as Jupiter notebooks, with application examples of the GLFM, including missing data estimation (a.k.a. matrix completion) and data exploratory analysis.

As an example, the script [`demo_toyImages'](https://github.com/ivaleraM/GLFM/blob/master/demos/python/demo_toy_images.ipynb) replicates the example of the IBP linear-Gaussian model in (Griffiths and Ghahramani, 2011) by generating a small set of images composed by different combinations of four original images plus additive Gaussian noise. Using the GLFM, we are able to recover the original images seamlessly.

Other examples include demo_matrix_completion_MNIST, demo_data_exploration_counties, and demo_data_exploration_prostate, available for PYTHON, Matlab and R. For more detail, please visit our [demo website](https://ivaleram.github.io/GLFM/demos.html). 

-------
Contact
-------

For further information or contact:

    Isabel Valera: isabel.valera.martinez (at) gmail.com
    Melanie F. Pradier: melanie.fpradier (at) gmail.com
    Maria Lomeni: mdl40  (at) cam.ac.uk


