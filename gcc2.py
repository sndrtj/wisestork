from math import *
import sys
import argparse
import statsmodels.nonparametric.smoothers_lowess as statlow
import numpy as np
from pyfaidx import Fasta

from countgc2 import get_gc_for_bin, get_n_per_bin
from utils import BedLine, Bin, attempt_integer


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


if __name__ == "__main__":
    desc = """
    Apply GC correction to sample BED files.
    """

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-I", "--input", type=str, required=True, help="Input BED file (produced by consam2.py)")
    parser.add_argument("-R", "--reference", type=str, required=True, help="Reference fasta (must be indexed)")
    parser.add_argument("-O", "--output", type=str, required=True, help="Path to output file")

    parser.add_argument("-n", "--frac-n", type=float, default=0.1, help="Maximum fraction of N-bases per bin")
    parser.add_argument("-r", "--frac-r", type=float, default=0.0001, help="Minimum fraction of reads per bin")
    parser.add_argument("-t", "--iter", type=int, default=3, help="Number of iterations for LOWESS function")
    parser.add_argument("-l", "--frac-lowess", type=float, default=0.1, help="Fraction of data to use for LOWESS function")

    args = parser.parse_args()
    fasta = Fasta(args.reference)
    bed_lines = [BedLine(*map(attempt_integer, x.split("\t"))) for x in open(args.input)]
    corrected = correct(bed_lines, fasta, args.frac_n, args.frac_r, args.iter, args.frac_lowess)

    with open(args.output, "wb") as ohandle:
        for line in corrected:
            ohandle.write(str(line) + "\n")



