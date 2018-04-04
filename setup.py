import os
import re
from setuptools import setup, find_packages
with open("deTiN/__about__.py") as reader:
    __version__ = re.search(
        r'__version__ ?= ?[\'\"]([\w.]+)[\'\"]',
        reader.read()
    ).group(1)


# Setup information
setup(
    name = 'deTiN',
    version = __version__,
    packages = find_packages(),
    description = 'Somatic analysis toolkit for dealing with tumor in normal contamination',
    author = 'Broad Institute - Cancer Genome Computational Analysis',
    author_email = 'amaro@broadinstitute.org',
    long_description = 'see publication',
    entry_points = {
        'console_scripts': [
            'deTiN = deTiN.deTiN:main'
        ]
    },
    install_requires = [
    'numpy',
    'matplotlib',
    'pandas',
    'scipy',
    'sklearn',
    'argparse'
    ],
    classifiers = [
        "Programming Language :: Python :: 2",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
)
