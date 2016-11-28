"""
wiseguy.count
~~~~~~~~~~~~~~

:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
:license: GPLv3
"""

import argparse
import pysam

from .utils import get_bins


def reads_per_bin(bam_reader, chromosome, bin):
    """
    Number of reads per bin
    :param bam_reader: an instance of pysam.AlignmentFile
    :param chromosome: chromosome name
    :param bin: Bin namedtuple
    :return: integer
    """
    try:
        reads = bam_reader.count(chromosome, bin.start, bin.end)
    except ValueError:
        reads = 0

    return reads


def get_chromosomes_from_header(header):
    """
    Get list of chromosome names from bam headers
    :param header: bam header
    :return: list of tuples of (chromosome names, chromosome_length)
    """
    return [(x['SN'], x['LN']) for x in header['SQ']]


def get_chromosomes_from_fasta(fa):
    """
    Get chromosome names from fasta file
    :param fa: instance of pyfaidx.Fasta
    :return: list of tuples of (chromosome names, chromosome_lenght)
    """
    # TODO
    pass


def count(input, output, binsize, reference):
    """
    Main function for counting reads per bin
    :param input: Path to input BAM
    :param output: path to output BED
    :param binsize: binsize
    :param reference: path to reference fasta
    """
    # TODO: check whether chromosome names in reference match those in bam file
    samfile = pysam.AlignmentFile(input)
    chromosomes = get_chromosomes_from_header(samfile.header)
    with open(output, "wb") as ohandle:
        for ch, ln in chromosomes:
            for bin in get_bins(ln, binsize):
                val = reads_per_bin(samfile, ch, bin)
                bed = "{0}\t{1}\t{2}\t{3}".format(ch, bin.start, bin.end, val)
                ohandle.write(bytes(bed + "\n"))


