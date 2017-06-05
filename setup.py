from setuptools import setup

setup(name='lprec_python',
      # random metadata. there's more you can supploy
      author='Viktor Nikitin',
      version='0.1',

      # this is necessary so that the swigged python file gets picked up
      py_modules=['lprec'],
      packages = ['lprecmods'])
   
