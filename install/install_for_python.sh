#!/bin/bash

if ! conda > /dev/null 2>&1; then
    echo "A requirement for the GLFM library is to previously install Anaconda at: https://www.anaconda.com/download/"
    exit 1
fi

conda install -y gsl
conda install -y -c anaconda cython=0.25.2
conda install -y -c pesoto cython_gsl=1.0.0
conda install -y -c anaconda cymem=1.31.2
cd ../src/Ccode/wrapper_python
python setup.py build_ext --inplace
