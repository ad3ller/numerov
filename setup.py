#!python
from setuptools import setup

setup(name='numerov',
      version='0.0.2',
      description='numerov wf solver',
      author='Adam Deller',
      author_email='a.deller@ucl.ac.uk',
      license='BSD 3-clause',
      packages=['numerov'],
      install_requires=[
          'numpy'
      ],
      zip_safe=False)