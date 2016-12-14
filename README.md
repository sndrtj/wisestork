![travis](https://travis-ci.org/sndrtj/wisecondor)
![coveralls](https://coveralls.io/github/sndrtj/wisecondor?branch=master)
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

A typical workflow starts with BAM files. Those BAM files _must_ be
sorted and indexed. 

The first step in a Wiseguy analysis is the `count` step. This 
generates read counts per bin, and writes this to a BED file. The 
command to do this, would look like the following:

`wiseguy count -I <input.bam> -R <fasta.fa> -O <out.bed> -B <binszise>`
 
The `-B` flag can be left out: Wiseguy defaults to a binsize of 50kb.
However, you will likely want a different binsize.

Once you have the count BED file, we have to correct for GC bias. The
command to do this is:

`wiseguy gc-correct -I <input.bed> -R <fasta.fa> -O <out.gc.bed> -B <binsize>`

For the next step, we need the result bgzipped and tabixed, so you'll 
have to execute `bgzip <out.gc.bed> && tabix -pbed <out.gc.bed.gz>`

The last step, the `zscore` step, calculates Z-scores for each bin.
It requires you to have generated a reference dictionary beforehand. 
The command to create z-scores again looks pretty similar to the 
earlier two:

`wiseguy zscore -I <input.bed.gz> -R <fasta.fa> -O <out.z.bed> -D <dictionary.bed.gz> -B <binsize>`

### Creating reference dictionaries

The above assumes you have already created a reference dictionary. 
If this is not the case, you will have to generate this file. 

To create the reference dictionary you will need a set of gc-corrected
BED files (from `wiseguy gc-correct`) of normal samples, and feed those
to `wiseguy newref`. The rewref command will then find the nearest
neighbours of every bin. Later on, in the zscore command, this
information is used to get a set of "reference bins" from the query
sample. 

Command to be used:

`wiseguy newref -I <input.gz.bed> -I <input2.gz.bed> [...] -O <out.ref.bed> -R <fasta.fa> -B <binsize>`
 
The output of this _must_ be sorted with bedtools, and then bgzipped
and tabixed. 

### Usage

```
Usage: wiseguy [OPTIONS] COMMAND [ARGS]...

  Discover CNVs from BAM files.

  A typical workflow first extracts regions from a BAM file
  The resulting BED tracks must then be GC-corrected.
  Using a reference track of region similarity,
  One can then calculate Z-scores for every region.

  The following sub-commands are supported:
   - count: count coverage per bin
   - gc-correct: GC-correct bins
   - zscore: calculate Z-scores
   - newref: Generate a new reference dictionary of bin similarities

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  count       Count coverages
  gc-correct  GC correct
  newref      Create new reference
  zscore      Calculte Z-scores
```

You can additional help by typing `wiseguy <command> --help`

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


## TODO

* Stouffer test, call and plot functions.
* In-memory option for `zscore`, that should be much faster.
* Some general speed ups. 
* Option to supply a BED file with bins, as opposed to generating them
from the reference sequence.