import pysam
import pyfaidx

from wiseguy.count import *


class TestFunctions:

    def test_sam_header(self):
        sam = pysam.AlignmentFile("test/data/test.bam")
        chroms = get_chromosomes_from_header(sam.header)
        assert len(chroms) == 1
        assert chroms[0][0] == "chrQ"
        assert chroms[0][1] == 500

    def test_fq_header(self):
        fa = pyfaidx.Fasta("test/data/chrQ.fasta")
        chroms = get_chromosomes_from_fasta(fa)
        assert len(chroms) == 1
        assert chroms[0][0] == "chrQ"
        assert chroms[0][1] == 500

    def test_count(self):
        sam = pysam.AlignmentFile("test/data/test.bam")
        bins = get_bins(500, 100)
        assert reads_per_bin(sam, "chrQ", bins[0]) == 40
        assert reads_per_bin(sam, "chrQ", bins[1]) == 52
        assert reads_per_bin(sam, "chrQ", bins[2]) == 34
        assert reads_per_bin(sam, "chrQ", bins[3]) == 42
        assert reads_per_bin(sam, "chrQ", bins[4]) == 26
        assert reads_per_bin(sam, "NotExist", bins[0]) == 0
