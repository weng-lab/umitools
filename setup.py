"""A toolset for handling sequencing data with unique molecular identifiers (UMIs)
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

setup(name="umitools",
      version='0.1.0',
      description='A toolset for handling sequencing data with unique molecular identifiers (UMIs)',
      url='https://github.com/weng-lab/umitools',
      author="Yu Fu",
      author_email="yfu@bu.edu",
      license="",
      classifiers=['Development Status :: 3 - Alpha',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Intended Audience :: Science/Research',
                   'Intended Audience :: Developers',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 2.6',
                   'Programming Language :: Python :: 2.7',
                   ],
      keywords='umi umitools unique molecular identifier',
      packages=find_packages(exclude=['contrib', 'docs', 'tests'])
)
