Wiseguy
=======
This is a complete re-implementation of the original Wisecondor program.
Its original purpose was to detect trisomies and smaller CNVs in 
maternal plasma samples using low-coverage WGS.
 
Wiseguy adds practical support for small bin sizes,
and is intended to be useful on regular WGS and Exome sequencing as well

For a full overview of differences with the original Wisecondor, 
see section Differences.

## Installation

The following system dependencies are required

* Python 3.4+

Furthermore, the following python packages are required:

* numpy
* matplotlib
* biopython
* statsmodels
* sklearn
* pysam
* pyfaidx
* click

It is recommended you use a virtualenv. 

To install wiseguy, create a virtualenv, install the python 
requirements using `pip install -r requirements.txt` and then run
`python setup.py develop`


## Input 

Wiseguy takes BAM files as input. These BAM files _must_ be indexed.
 
Additionally, you must provide a reference Fasta file, which should
likewise be indexed with `samtools faidx <fasta>`.  

## Running

_TODO_

## Differences

There are several important differences between this re-implementation
and the original wisecondor. 
 
* This re-implementation is organized as a regular python package, 
  while exposing several command-line tools. 
* Python 3 support
* All command-line tools now have UNIX-style argument parsing
* Generating reference sets for small bin sizes is now possible in 
  much less time. 
* Pickle files are no longer used. The output format is now regular BED,
  with a possible additional column. This means results can be used by 
  common downstream tools like Bedtools.
* The countgc step is now redundant. Its functionality is now integrated
  in the gcc step. 
* The reference bin selection method was modified. The
  original wisecondor calculated differences for every bin against every
  bin of every sample, and then repeated this calculation for every 
  chromosome. As this is an exponential operation, this made 
  reference bin selection prohibitively slow and memory-consuming for 
  smaller bin size. In stead of calculating differences, the new method
  applies a method (e.g. median) over the same bins of all samples, 
  and then _sorts_ the resulting list of bins. Similar bins can be
  selected using regular list slicing. This means the time complexity 
  of creating a new reference set is now just loglinear. Additional 
  filterings were left the same. 
* Use of the `statsmodels` lowess function, rather than biopython's. 
  This results in a significant speed-up of the gc correction.
