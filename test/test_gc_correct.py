from pyfaidx import Fasta
from pytest import fixture

from wiseguy.gc_correct import filter_bin, correct
from wiseguy.utils import BedLine, attempt_numeric

import sys

@fixture
def bedlines():
    with open("test/data/count.bed") as count_handle:
        bedlines = [BedLine(*map(attempt_numeric, x.split("\t"))) for x in count_handle]
    return bedlines

@fixture
def fasta():
    return Fasta("test/data/chrQ.fasta")


class TestFunctions:
    def test_filter_bin(self, bedlines, fasta):
        assert all([filter_bin(x, fasta, 0.1, 0.001) for x in bedlines])
        assert all([not filter_bin(x, fasta, 0.1, 5) for x in bedlines])

    def test_correct(self, bedlines, fasta):
        assert all([x.value == 0 for x in correct(bedlines, fasta, frac_r=5)])
        corrected = correct(bedlines, fasta)
        for i in range(4):
            assert 0.9 < corrected[i].value < 1.1
        assert 0.4 < corrected[4].value < 0.6



