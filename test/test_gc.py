from pyfaidx import Fasta

from wiseguy.gc import *
from wiseguy.utils import get_bins


class TestFunctions:

    def test_n_per_bin(self):
        fa = Fasta("test/data/chrQ.fasta")
        bins = get_bins(500, 100)
        for bin in bins:
            assert get_n_per_bin(fa, 'chrQ', bin) == 0

    def test_get_gc_for_bin(self):
        fa = Fasta("test/data/chrQ.fasta")
        bins = get_bins(500, 100)
        assert get_gc_for_bin(fa, "chrQ", bins[0]) == 47
        assert get_gc_for_bin(fa, "chrQ", bins[1]) == 50
        assert get_gc_for_bin(fa, "chrQ", bins[2]) == 54
        assert get_gc_for_bin(fa, "chrQ", bins[3]) == 58
        assert get_gc_for_bin(fa, "chrQ", bins[4]) == 50