#!python
from setuptools import setup, Extension
import numpy

try:
    from Cython.Build import cythonize
    ext_modules = cythonize(Extension("core", ["numerov/cy/core.pyx"], include_dirs=[numpy.get_include()]),
                            compiler_directives={'language_level' : 3, 'embedsignature': True})
except ImportError:
    ext_modules = []
except:
    raise


setup(name='numerov',
      version='0.0.5',
      description='numerov wf',
      author='Adam Deller',
      author_email='a.deller@ucl.ac.uk',
      license='BSD 3-clause',
      packages=['numerov'],
      install_requires=[
          'numpy'
      ],
      ext_package='numerov/cy',
      ext_modules=ext_modules,
      zip_safe=False)
