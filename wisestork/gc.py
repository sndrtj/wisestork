"""
wisestork.gc
~~~~~~~~~~

:copyright: (c) 2013 Roy Straver
:copyright: (c) 2013 VU University Medical Center
:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
:license: GPLv3
"""

from Bio.SeqUtils import GC


def get_gc_for_bin(fasta, chromosome, bin):
    """
    Get number of GC bases in a region
    :param fasta: an instance of pyfaidx.Fasta
    :param chromosome: chromosome name
    :param bin: Bin namedtuple
    :return: integer
    """
    perc = GC(fasta[chromosome][bin.start:bin.end].seq)
    dist = bin.end-bin.start
    return int((perc*dist)/100)


def get_n_per_bin(fasta, chromosome, bin):
    """
    Get number of N bases in a region
    :param fasta: an instance of pyfaidx.Fasta
    :param chromosome: chromosome name
    :param bin: Bin namedtuple
    :return: integer
    """
    return fasta[chromosome][bin.start:bin.end].seq.upper().count('N')
