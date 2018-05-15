"""A toolset for handling sequencing data with unique molecular identifiers (UMIs)
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
from os import path

setup(name="umitools",
      # Suffix a-z for PyPI test
      # Remove the suffix before uploading to PyPI
      version='0.3.3',
      description='A toolset for handling sequencing data with unique '
      'molecular identifiers (UMIs)',
      # Now explicitly requiring Python 3
      python_requires='>=3',
      url='https://github.com/weng-lab/umitools',
      author="Yu Fu",
      author_email="yfu@bu.edu",
      license="GPLv3",
      classifiers=['Development Status :: 3 - Alpha',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Intended Audience :: Science/Research',
                   'Intended Audience :: Developers',
                   'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
                   # 'Programming Language :: Python :: 2',
                   # 'Programming Language :: Python :: 2.6',
                   # 'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3',
                   'Programming Language :: Python :: 3.4',
                   'Programming Language :: Python :: 3.5',
                   'Programming Language :: Python :: 3.6',
                   ],
      # download_url='https://github.com/weng-lab/umitools/archive/0.2.0.tar.gz',
      keywords=" ".join(["RNA-seq", "small RNA-seq",
                         "unique molecular identifier",
                         "UMI", "PCR duplicates", "PCR cycle",
                         "starting material", "sequencing depth",
                         "high-throughput sequencing",
                         "deep sequencing", "transcriptome",
                         "genome", "umitools", "RNA", "small RNA",
                         "sequencing"]),
      install_requires=['pysam', 'editdistance', 'networkx'],
      packages=find_packages(exclude=['contrib', 'docs', 'tests']),
      # packages=['umitools'],
      # data_files=[('my_data', ['data/umitools.test.r1.fq.gz', 'data/umitools.test.r2.fq.gz'])],
      py_modules=["umitools/umi"],
      package_data={
          "sample": ['data/*.fq.gz'],
          },
      entry_points={
          'console_scripts': [
              'umitools=umitools.umitools:main',
              'umi_reformat_fastq=umitools.reformat_umi_fastq:main',
              'umi_reformat_sra_fastq=umitools.reformat_umi_sra_fastq:main',
              'umi_mark_duplicates=umitools.umi_mark_duplicates:main',
              'umi_simulator=umitools.umi_simulator:main'
              ]
              }
      # scripts=['scripts/reformat_umi_fastq', 'scripts/umi_mark_duplicates'],
)


