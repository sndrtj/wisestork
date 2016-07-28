import argparse
from Bio.SeqUtils import GC
from pyfaidx import Fasta

from utils import get_bins  # TODO : make package so we can do relative import


def get_gc_for_bin(fasta, chromosome, bin):
    """
    Get number of GC bases in a region
    :param fasta: an instance of pyfaidx.Fasta
    :param chromosome: chromosome name
    :param bin: 2-tuple of (start, end)
    :return: integer
    """
    perc = GC(fasta[chromosome][bin[0]:bin[1]].seq)
    return int((bin[1]-bin[0])*(perc/100))


if __name__ == "__main__":
    desc = """
    Count binned GC content in a fasta file.
    Outputs bed file with 4 column as number of GC-bases per bin.
    Your fasta file must be in indexed with e.g. samtools
    """
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-I", "--input", required=True, type=str, help="Input fasta file")
    parser.add_argument("-O", "--output", required=True, type=str, help="Output bed file")
    parser.add_argument("-b", "--binsize", type=int, default=int(1e6), help="Binsize")

    args = parser.parse_args()

    fasta = Fasta(args.input)
    with open(args.output, "w") as ohandle:
        for chromosome in fasta.keys():
            for bin in get_bins(len(fasta[chromosome]), args.binsize):
                val = get_gc_for_bin(fasta, chromosome, bin)
                bed = "{0}\t{1}\t{2}\t{3}".format(chromosome, bin[0], bin[1], val)
                ohandle.write(bed + "\n")