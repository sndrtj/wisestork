"""
wisestork.gc_correct

:copyright: (c) 2013 Roy Straver
:copyright: (c) VU University Medical Center
:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
:license: GPLv3
"""

import numpy as np
import warnings

import statsmodels.nonparametric.smoothers_lowess as statlow
from pyfaidx import Fasta

from .utils import BedLine, attempt_numeric
from .gc import get_gc_for_bin, get_n_per_bin


def filter_bin(bin, ref_fasta, frac_n, frac_r):
    """
    Return True if a 'correct' bin, return False if an 'incorrect' bin
    :param bin: BedLine namedtuple
    :param ref_fasta: an instance of pyfaidx.Fasta
    :param frac_n: maximal fraction of N-bases per bin
    :param frac_r: minimum fraction of reads per bin
    :return: Boolean
    """
    ns = get_n_per_bin(ref_fasta, bin.chromosome, bin)
    return ns < ((bin.end - bin.start)*frac_n) and bin.value > ((bin.end - bin.start)*frac_r)  # noqa


def correct(inputs, fasta, frac_n=0.1, frac_r=0.0001, lowess_iter=3,
            lowess_frac=0.1):
    """
    GC-correct input bed lines.
    GC correction takes place with a local regression (LOWESS) on GC perc vs
    number of reads
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
        if filter_bin(line, fasta, frac_n, frac_r):
            gcs.append(get_gc_for_bin(fasta, line.chromosome, line))
            reads.append(line.value)

    reads = np.array(reads, np.float)
    gcs = np.array(gcs, np.float)
    if lowess_frac*len(reads) < 4 and len(reads) > 0:
        # need at least four data ponts
        warnings.warn("Too few data points for lowess. Raising lowess_frac")
        lowess_frac = 4.0/len(reads)
        delta = 0  # remove delta in this case
    else:
        delta = 0.01 * len(gcs)
    lowess = statlow.lowess(reads, gcs, return_sorted=False,
                            delta=delta, frac=lowess_frac,
                            it=lowess_iter).tolist()

    corrected_lines = []

    for line in inputs:
        if filter_bin(line, fasta, frac_n, frac_r):
            corr_val = float(line.value) / lowess.pop(0)
        else:
            corr_val = 0
        n_bed = BedLine(line.chromosome, line.start, line.end, corr_val)
        corrected_lines.append(n_bed)

    return corrected_lines


def gc_correct(input, output, reference, frac_n, frac_r, iter, frac_lowess):
    fasta = Fasta(reference)
    bed_lines = [BedLine(*map(attempt_numeric, x.split("\t"))) for
                 x in open(input)]
    corrected = correct(bed_lines, fasta, frac_n, frac_r, iter, frac_lowess)

    with open(output, "wb") as ohandle:
        for line in corrected:
            ohandle.write(bytes(str(line) + "\n", 'utf-8'))
