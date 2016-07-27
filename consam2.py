import argparse
import pysam


def get_bins(chromosome_length, binsize):
    """
    Get list of 2-tuples of start and end positions of bins for a given chromosome
    :param chromosome_length: integer
    :param binsize: integer
    :return: list of 2-tuples of (start, end). Start = 0-based
    """

    starts = list(range(0, chromosome_length, binsize))
    ends = []
    for s in starts:
        if s + binsize < chromosome_length:
            ends.append(s+binsize)
        else:
            ends.append(chromosome_length)
    return [(s, e) for s, e in zip(starts, ends)]


def bin_to_bed(bam_reader, chromosome, bin):
    """
    Create bed line from bin
    :param bam_reader: an instance of pysam.AlignmentFile
    :param chromosome: chromosome name
    :param bin: 2-tuple of (start, end)
    :return: bedline of chr\tstart\tend\tnreads
    """
    try:
        reads = bam_reader.count(chromosome, bin[0], bin[1])
    except ValueError:
        reads = 0

    return "{0}\t{1}\t{2}\t{3}".format(chromosome, bin[0], bin[1], reads)


def get_chromosomes_from_header(header):
    """
    Get list of chromosome names from bam headers
    :param header: bam header
    :return: list of tuples of (chromosome names, chromosome_length)
    """
    return [(x['SN'], x['LN']) for x in header['SQ']]


if __name__ == "__main__":
    desc = """
    This script takes a BAM file, and calculates the number of reads per bin.
    It will output a BED file (with 0-based positions) of regions and associated reads/bin.

    Your BAM file must be indexed with tabix, and must contain chromosome names and lengths in the header.
    """
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-I", "--input", type=str, required=True, help="Input BAM file")
    parser.add_argument("-O", "--output", type=str, required=True, help="Output BED file")

    parser.add_argument("-b", "--binsize", type=int, default=int(1e6), help="Binsize")

    args = parser.parse_args()

    samfile = pysam.AlignmentFile(args.input)
    chromosomes = get_chromosomes_from_header(samfile.header)
    with open(args.output, "wb") as ohandle:
        for ch, ln in chromosomes:
            for bin in get_bins(ln, args.binsize):
                bed = bin_to_bed(samfile, ch, bin)
                ohandle.write(bed + "\n")


