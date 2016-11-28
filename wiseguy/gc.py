"""
wiseguy.gc
~~~~~~~~~~

:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
:license: GPLv3
"""

import argparse
from Bio.SeqUtils import GC
from pyfaidx import Fasta

from .utils import get_bins


def get_gc_for_bin(fasta, chromosome, bin):
    """
    Get number of GC bases in a region
    :param fasta: an instance of pyfaidx.Fasta
    :param chromosome: chromosome name
    :param bin: Bin namedtuple
    :return: integer
    """
    perc = GC(fasta[chromosome][bin.start:bin.end].seq)
    return int((bin.start-bin.end)*(perc/100))


def get_n_per_bin(fasta, chromosome, bin):
    """
    Get number of N bases in a region
    :param fasta: an instance of pyfaidx.Fasta
    :param chromosome: chromosome name
    :param bin: Bin namedtuple
    :return: integer
    """
    return fasta[chromosome][bin.start:bin.end].seq.upper().count('N')


if __name__ == "__main__":
    desc = """
    Count binned GC content in a fasta file.
    Outputs bed file with 4 column as number of GC-bases per bin.
    It will output number of hard-masked (= N) bases per bin in separate file
    Your fasta file must be in indexed with e.g. samtools
    """
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-I", "--input", required=True, type=str, help="Input fasta file")
    parser.add_argument("-O", "--output", required=True, type=str, help="Output bed file")
    parser.add_argument("-M", "--hard-masked", required=True, type=str, help="Output hard-masked file")
    parser.add_argument("-b", "--binsize", type=int, default=int(1e6), help="Binsize")

    args = parser.parse_args()

    fasta = Fasta(args.input)
    with open(args.output, "w") as ohandle, open(args.hard_masked, "wb") as masked_handle:
        for chromosome in fasta.keys():
            for bin in get_bins(len(fasta[chromosome]), args.binsize):
                val = get_gc_for_bin(fasta, chromosome, bin)
                bed = "{0}\t{1}\t{2}\t{3}".format(chromosome, bin.start, bin.end, val)
                ohandle.write(bed + "\n")
                n_val = get_n_per_bin(fasta, chromosome, bin)
                m_bed= "{0}\t{1}\t{2}\t{3}".format(chromosome, bin.start, bin.end, n_val)
                masked_handle.write(m_bed + "\n")