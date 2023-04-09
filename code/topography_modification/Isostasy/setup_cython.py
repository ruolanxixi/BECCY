# Description: Setup file for Cython functions
#
# Building: python setup_cython.py build_ext --inplace
#
# Copyright (c) 2022 ETH Zurich, Christian R. Steger

# Load modules
from distutils.core import setup
from Cython.Distutils import build_ext
from distutils.extension import Extension
import numpy

ext_modules = [Extension("isostasy_cy",
               ["isostasy_cy.pyx"],
               libraries=["m", "iomp5", "pthread"],
               extra_compile_args=["-O3", "-ffast-math", "-fopenmp"],
               extra_link_args=["-fopenmp"],
               include_dirs=[numpy.get_include()])]

setup(name="isostasy_cy",
      cmdclass={"build_ext": build_ext},
      ext_modules=ext_modules)
