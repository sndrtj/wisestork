"""
setup.py
~~~~~~~~~~~~~
:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
:license: GPL
"""
from os.path import abspath, dirname, join

from setuptools import setup, find_packages


readme_file = join(abspath(dirname(__file__)), "README.md")
with open(readme_file) as desc_handle:
    long_desc = desc_handle.read()


setup(
    name="wisestork",
    version="0.1.0",
    description="",
    author="Sander Bollen",
    license="GPL",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "matplotlib",
        "statsmodels",
        "pysam",
        "pyfaidx",
        "biopython",
        "scipy",
        "click",
        "progressbar2"
    ],
    entry_points={
        "console_scripts": [
            "wisestork = wisestork.wisestork:main"
        ]
    },
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)

