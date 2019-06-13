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
    description="Within-sample CNV calling",
    author="Sander Bollen",
    author_email="sander@sndrtj.eu",
    python_requires=">=3.5",
    license="GPLv3+",
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
        "License :: OSI Approved :: GNU General Public License v3 or "
        "later (GPLv3+)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ]
)
