from distutils.core import setup, Extension
import os
import numpy

#python_include_dir = '/home/joao/anaconda/include'
numpy_include_dir = os.path.dirname(numpy.__file__) + '/core/include/numpy'
# numpy_include_dir = '/home/joao/anaconda/lib/python2.7/site-packages/numpy/core/include/numpy'

module1 = Extension('periodogram_CLEAN', 
        include_dirs = [numpy_include_dir],
        libraries = ['m'],
        sources = ['periodogram_CLEAN.c'])

setup (name = 'CLEAN algorithm',
        version = '1.0',
        description = '',
        ext_modules = [module1])
