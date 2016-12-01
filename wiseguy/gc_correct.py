"""
wiseguy.gc_correct

:copyright: (c) 2013 Roy Straver
:copyright: (c) VU University Medical Center
:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden Univeristy Medical Center
:license: GPLv3
"""

import argparse
import numpy as np

import statsmodels.nonparametric.smoothers_lowess as statlow
from pyfaidx import Fasta

from .utils import BedLine, Bin, attempt_integer
from .gc import get_gc_for_bin, get_n_per_bin


def filter_bin(value, chromosome, bin, ref_fasta, frac_n, frac_r):
    """
    Return True if a 'correct' bin, return False if an 'incorrect' bin
    :param value: number of reads in this bin
    :param chromosome: name of chromosome
    :param bin: Bin namedtuple
    :param ref_fasta: an instance of pyfaidx.Fasta
    :param frac_n: maximal fraction of N-bases per bin
    :param frac_r: minimum fraction of reads per bin
    :return: Boolean
    """
    ns = get_n_per_bin(ref_fasta, chromosome, bin)
    return ns < ((bin.end - bin.start)*frac_n) and value > ((bin.end - bin.start)*frac_r)


def correct(inputs, fasta, frac_n=0.1, frac_r=0.0001, lowess_iter=3, lowess_frac=0.1):
    """
    GC-correct input bed lines.
    GC correction takes place with a local regression (LOWESS) on GC perc vs number of reads
    :param inputs: list of BedLine namedtuples
    :param fasta: instance of pyfaidx.Fasta
    :param frac_n: maximal fraction on N-bases per bin
    :param frac_r: minimum fraction of reads per bin
    :param lowess_iter: amount of iterations of LOWESS function
    :param lowess_frac: fraction of input data used for LOWESS function
    :return: corrected BedLines
    """
    reads = []
    gcs = []
    for line in inputs:
        chromosome = line.chromosome
        bin = Bin(line.start, line.end)
        if filter_bin(line.value, chromosome, bin, fasta, frac_n, frac_r):
            gcs.append(get_gc_for_bin(fasta, chromosome, bin))
            reads.append(line.value)

    reads = np.array(reads, np.float)
    gcs = np.array(gcs, np.float)
    delta = 0.01 * len(gcs)
    lowess = statlow.lowess(reads, gcs, return_sorted=False,
                            delta=delta, frac=lowess_frac,
                            it=lowess_iter).tolist()

    corrected_lines = []

    for line in inputs:
        chromosome = line.chromosome
        bin = Bin(line.start, line.end)
        if filter_bin(line.value, chromosome, bin, fasta, frac_n, frac_r):
            corr_val = float(line.value) / lowess.pop(0)
        else:
            corr_val = 0
        n_bed = BedLine(chromosome, line.start, line.end, corr_val)
        corrected_lines.append(n_bed)

    return corrected_lines


def gc_correct(input, output, reference, frac_n, frac_r, iter, frac_lowess):
    fasta = Fasta(reference)
    bed_lines = [BedLine(*map(attempt_integer, x.split("\t"))) for x in open(input)]
    corrected = correct(bed_lines, fasta, frac_n, frac_r, iter, frac_lowess)

    with open(output, "wb") as ohandle:
        for line in corrected:
            ohandle.write(bytes(str(line) + "\n", 'utf-8'))



