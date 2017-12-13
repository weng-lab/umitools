"""A toolset for handling sequencing data with unique molecular identifiers (UMIs)
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

setup(name="umitools",
      version='0.2.0',
      description='A toolset for handling sequencing data with unique molecular identifiers (UMIs)',
      url='https://github.com/weng-lab/umitools',
      author="Yu Fu",
      author_email="yfu@bu.edu",
      license="GPLv3",
      classifiers=['Development Status :: 3 - Alpha',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Intended Audience :: Science/Research',
                   'Intended Audience :: Developers',
                   'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 2.6',
                   'Programming Language :: Python :: 2.7',
                   ],
      download_url='https://github.com/weng-lab/umitools/archive/0.2.0.tar.gz', 
      keywords=['umi', 'umitools', 'unique molecular identifier', 'RNA', 'small RNA', 'sequencing'],
      install_requires=['pysam', 'editdistance', 'networkx'],
      packages=find_packages(exclude=['contrib', 'docs', 'tests']),
      # data_files=[('my_data', ['data/umitools.test.r1.fq.gz', 'data/umitools.test.r2.fq.gz'])],
      package_data={
          "umitools": ['data/*.fq.gz'],
          },
      entry_points={
          'console_scripts': [
              'reformat_fastq=umitools.reformat_umi_fastq:main',
              'mark_duplicates=umitools.umi_mark_duplicates:main',              
              ]
              }
      # scripts=['scripts/reformat_umi_fastq', 'scripts/umi_mark_duplicates'],
)
