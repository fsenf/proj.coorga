

from setuptools import setup, find_packages

setup(name='coorga',
      version='0.1',
      description='Collection of Convective Organization Tools',
      author='Fabian Senf',
      author_email='senf@tropos.de',
      license='GPL',
      packages=find_packages(),
#['coorga', 'coorga.inout', 'coorga.object_creation', 'coorga.metrics'],
      zip_safe=False)

