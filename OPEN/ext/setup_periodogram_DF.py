from distutils.core import setup, Extension
import os
import numpy

#python_include_dir = '/home/joao/anaconda/include'
numpy_include_dir = os.path.dirname(numpy.__file__) + '/core/include/numpy'
# numpy_include_dir = '/home/joao/anaconda/lib/python2.7/site-packages/numpy/core/include/numpy'

module1 = Extension('periodogram_DF', 
        include_dirs = [numpy_include_dir],
        libraries = ['m'],
        sources = ['periodogram_DF.c'])

setup (name = 'Lomb-Scargle periodogram with phases output',
        version = '1.0',
        description = '',
        ext_modules = [module1])
