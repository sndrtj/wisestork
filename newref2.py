import argparse
import numpy as np

from utils import BedLine, get_bins
from pyfaidx import Fasta


def get_unique_bins(fasta, binsize):
    """
    Get a list of unique bins (by position, not value!)
    :param fasta: pyfaidx.Fasta instance
    :return: bins
    """
    bins = []
    for contig in fasta:
        n = contig.name
        tmp_bins = get_bins(len(contig), binsize)
        bins += [BedLine(n, x[0], x[1], 0) for x in tmp_bins]

    return bins


def build_main_list(bins, binsize, reference):
    """
    Build main list of bins.
    All unique positions are selected
    Then the median value of all bins on those positions are calculated
    :param bins: list of bins
    :param binsize: binsize:
    :param reference: instance of pyfaidx.Fasta
    :return: list of bins (1 per position), sorted by median value
    """
    unique_positions = get_unique_bins(reference, binsize)
    median_bins = []
    for i, u in enumerate(unique_positions):
        working_bins = bins[i::len(unique_positions)]
        med = np.median([x.value for x in working_bins])
        median_bins.append(BedLine(u.chromosome, u.start, u.end, med))

    return sorted(median_bins, key=lambda x: x.value)


class ReferenceBinGenerator(object):
    """
    Iterator for generating reference bins
    Every item will be 2-tuple of (Bin, List[Bins])
    With the second item being similar bins
    """

    def __init__(self, inputs, n_bins, reference, binsize=int(1e6)):
        """
        Create instance of ReferenceBinGenerator
        :param inputs: list of paths to files of gc-corrected bedgraph files
        :param n_bins: number of neighbour bins to consider
        """
        self.inputs = inputs
        self.n_bins = n_bins
        self.fasta = Fasta(reference)
        self.binsize = binsize
        self.__bins = self.get_all_bins()
        self.__idx = 0

    def get_all_bins(self):
        bins = []
        for inp in self.inputs:
            with open(inp) as ihandle:
                for record in ihandle:
                    vals = record.strip().split("\t")
                    bins += [BedLine(vals[0], int(vals[1]), int(vals[2]), float(vals[3]))]
        return build_main_list(bins, self.binsize, self.fasta)

    def __next__(self):
        if self.__idx == len(self.__bins):
            raise StopIteration
        nearest = self.get_nearest_at_idx(self.__idx)
        cur = self.__bins[self.__idx]
        filtered = self.filter_bins(cur, nearest)
        self.__idx += 1
        return cur, filtered

    def next(self):
        return self.__next__()

    def __iter__(self):
        return self

    def get_nearest_at_idx(self, idx):
        if idx < self.n_bins:
            vals = self.__bins[:idx]
        elif idx > len(self.__bins) - self.n_bins:
            vals = self.__bins[idx:]
        else:
            vals = self.__bins[idx-(self.n_bins/2):idx+(self.n_bins/2)]
        return vals

    def filter_bins(self, target_bin, bins):
        """
        Filter a list of reference bins
        Bins to exclude:
         * target bin
         * Bins with value > 3x mean+stdev
         * Adjacent bins # TODO
        :param target_bin:
        :param bins:
        :return: filtered list of bins (may be empty)
        """
        ref_bins = [x for x in bins if
                    x.chromosome != target_bin.chromosome and
                    x.start != target_bin.start and
                    x.end != target_bin.end]
        if len(ref_bins) == 0:
            return ref_bins

        stdev = np.std([x.value for x in ref_bins])
        mean = np.mean([x.value for x in ref_bins])
        non_outliers = [x for x in ref_bins if mean+(3*stdev) > x.value > mean-(3*stdev)]
        return non_outliers


if __name__ == "__main__":
    desc = """
    Create a reference database from a set of GC-corrected bedgraph files.

    Your bedgraph file(s) _must_ be sorted in the order of the reference fasta.
    You should sort the resulting file with bedtools, bgzip and tabix it.
    E.g.
    python newref2.py <<arguments>> -O out.bed && bedtools sort -i out.bed | bgzip -c > out.bed.gz && tabix -pbed out.bed.gz
    """

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-I", "--input", help="Input file(s) for reference building", action="append", required=True)
    parser.add_argument("-O", "--output", help="Output bed-likefile", required=True)
    parser.add_argument("-R", "--reference", help="Reference fasta (must be indexed)", required=True)
    parser.add_argument("-b", "--binsize", help="binsize", type=int, default=int(1e6))
    parser.add_argument("--n-bins", help="Number of reference bins to consider for each bin", type=int, default=250)

    args = parser.parse_args()

    gen = ReferenceBinGenerator(args.input, args.n_bins, args.reference, args.binsize)
    refs = []
    ohandle = open(args.output, "w")
    for x in gen:
        chrom = x[0].chromosome
        start = x[0].start
        end = x[0].end
        references = "|".join([str(x).replace("\t", ",") for x in x[1]])
        if len(references) == 0:
            references = np.nan

        t = BedLine(chrom, start, end, references)
        ohandle.write(str(t) + "\n")
