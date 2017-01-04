"""
setup.py
~~~~~~~~~~~~~
:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
:license: GPL
"""
from os.path import abspath, dirname, join

from setuptools import setup

from wiseguy import __version__


readme_file = join(abspath(dirname(__file__)), "README.md")
with open(readme_file) as desc_handle:
    long_desc = desc_handle.read()


setup(
    name="wiseguy",
    version=__version__,
    description="",
    author="Sander Bollen",
    license="GPL",
    install_requires=[
        "numpy",
        "matplotlib",
        "statsmodels",
        "pysam",
        "pyfaidx",
        "biopython",
        "click"
    ],
    entry_points={
        "console_scripts": [
            "wiseguy = wiseguy.wiseguy:main"
        ]
    },
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)

