[![Build Status](https://travis-ci.org/sndrtj/wisecondor.svg?branch=master)](https://travis-ci.org/sndrtj/wisecondor)
[![codecov](https://codecov.io/gh/sndrtj/wisecondor/branch/master/graph/badge.svg)](https://codecov.io/gh/sndrtj/wisecondor)

Wisestork
=======
This is a complete re-implementation of the original 
[Wisecondor](https://github.com/VUmcCGP/wisecondor) program.
Its original purpose was to detect trisomies and smaller CNVs in 
maternal plasma samples using low-coverage WGS.
 
Wisestork adds practical support for small bin sizes,
and is intended to be useful on regular WGS and Exome sequencing as well.

For a full overview of differences with the original Wisecondor, 
see section Differences.

## Installation

The following system dependencies are required

* Python 3.5+

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

To install wisestork, create a virtualenv, install the python 
requirements using `pip install -r requirements.txt` and then run
`python setup.py develop`


## Input 

Wisestork takes BAM files as input. These BAM files _must_ be indexed.
 
Additionally, you must provide a reference Fasta file, which should
likewise be indexed with `samtools faidx <fasta>`.  

## Running

A typical workflow starts with BAM files. Those BAM files _must_ be
sorted and indexed. 

The first step in a Wisestork analysis is the `count` step. This 
generates read counts per bin, and writes this to a BED file. The 
command to do this, would look like the following:

`wisestork count -I <input.bam> -R <fasta.fa> -O <out.bed> -B <binszise>`
 
The `-B` flag can be left out: Wisestork defaults to a binsize of 50kb.
However, you will likely want a different binsize.

Once you have the count BED file, we have to correct for GC bias. The
command to do this is:

`wisestork gc-correct -I <input.bed> -R <fasta.fa> -O <out.gc.bed> -B <binsize>`

For the next step, we need the result bgzipped and tabixed, so you'll 
have to execute `bgzip <out.gc.bed> && tabix -pbed <out.gc.bed.gz>`

The last step, the `zscore` step, calculates Z-scores for each bin.
It requires you to have generated a reference dictionary beforehand. 
The command to create z-scores again looks pretty similar to the 
earlier two:

`wisestork zscore -I <input.bed.gz> -R <fasta.fa> -O <out.z.bed> -D <dictionary.bed.gz> -B <binsize>`


### User-supplied bins

In stead of supplying a bin _size_ for each step, you may also supply a 
bin _file_. This file should be a (preferably sorted) BED file with regions
that exist in the input BAM file. This option is primarily useful for 
WES analyses, where the bin file would correspond to a target/bait region
file. Please do note that contigs must be identical to those in the 
input BAM file. 

You can supply a bin file using the `-L` flag for any subcommand.
This will supersede any usage of the `-B` flag.

### Creating reference dictionaries

The above assumes you have already created a reference dictionary. 
If this is not the case, you will have to generate this file. 

To create the reference dictionary you will need a set of gc-corrected
BED files (from `wisestork gc-correct`) of normal samples, and feed those
to `wisestork newref`. The rewref command will then find the nearest
neighbours of every bin. Later on, in the zscore command, this
information is used to get a set of "reference bins" from the query
sample. 

Command to be used:

`wisestork newref -I <input.gz.bed> -I <input2.gz.bed> [...] -O <out.ref.bed> -R <fasta.fa> -B <binsize>`
 
The output of this _must_ be sorted with bedtools, and then bgzipped
and tabixed. 

### Usage

```
Usage: wisestork [OPTIONS] COMMAND [ARGS]...

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
  zscore      Calculate Z-scores
```

You can additional help by typing `wisestork <command> --help`

## Differences

There are several important differences between this re-implementation
and the original wisecondor. 
 
* This re-implementation is organized as a regular python package, 
  while exposing several command-line tools. 
* Python 3 support. In fact, it's only tested on python 3.
* All command-line tools now have UNIX-style argument parsing
* Generating reference sets for small bin sizes is now possible in 
  much less time. 
* Pickle files are no longer used. The output format is now regular BED,
  with a possible additional column. This means results can be used by 
  common downstream tools like Bedtools.
* User supplied bin files in regular BED format. 
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

## Naming

Why name this tool wisestork, you might think?
Well, a _condor_ is a bird. As this is a re-implementation / fork of
wisecondor, I figured another bird would be nice name. As I live in The Hague,
and The Hague has a stork as a city symbol, I put one and one together.
Thus, _wisestork_ was born. 

## License

GPLv3