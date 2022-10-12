#!/usr/bin/env python
import sys
from setuptools import setup, Extension
from distutils.errors import DistutilsError
from setuptools.command.build_py import build_py as _build
from setuptools.command.develop import develop as _develop
import subprocess as sp
import os, sys


# Get numpy include directory (works across versions)
import numpy
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

libs = ['cute']#'sharp2', 'cfitsio', 'gsl', 'gslcblas', 'm'] + FFTW_LIBS

_cutelib = Extension("_cutelib",
                     ["cutepy/cute_wrap.c"],
                     extra_objects=[],
                     libraries=libs,
                     #library_dirs=["./_deps/lib/"],
                     include_dirs=[numpy_include, "./src/"],#,"./_deps/include/"],
                     #extra_compile_args=extra,
                     #extra_link_args=extra
)

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name="cutepy",
      version="0.0",
      author="David Alonso",
      author_email="david.alonso@physics.ox.ac.uk",
      description="Library for angular correlations",
      long_description=long_description,
      long_description_content_type="text/markdown",
      url="https://github.com/CutePy",
      #cmdclass={'build_py': build, 'develop': develop},
      classifiers=[
          'License :: OSI Approved :: BSD License',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.6',
          'Operating System :: Unix',
          'Operating System :: MacOS'],
      packages=['cutepy'],
      ext_modules=[_cutelib],
      )
