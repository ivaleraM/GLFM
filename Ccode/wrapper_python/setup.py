from distutils.core import setup
from Cython.Distutils import Extension
from Cython.Distutils import build_ext
import cython_gsl
import numpy

# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
import distutils.sysconfig
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

setup(
    include_dirs = [cython_gsl.get_include(), numpy.get_include()],
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("gsl_run",
            sources=["gsl_run.pyx","../core/GeneralFunctions.cpp","../core/InferenceFunctions.cpp"],
            language="c++",
            libraries=cython_gsl.get_libraries(),
            library_dirs=[cython_gsl.get_library_dir()],
            include_dirs=[cython_gsl.get_cython_include_dir()])]
    )
