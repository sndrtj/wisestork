import pytest
from wiseguy.utils import *

@pytest.fixture
def test_bedline():
    return BedLine(b"chr1", 1, 100, 20)

class TestFunctions:

    def test_utf8(self):
        assert utf8(b"test") == b"test"
        assert utf8("test") == b"test"
        assert utf8(1) == b"1"

    def test_as_str(self):
        assert as_str("test") == "test"
        assert as_str(b"test") == "test"

    def test_attempt_numeric(self):
        assert attempt_numeric("2") == 2
        assert attempt_numeric(b"2") == 2
        assert attempt_numeric("2.5") == 2.5
        assert attempt_numeric(b"2.5") == 2.5
        assert attempt_numeric("test") == "test"
        assert attempt_numeric(b"test") == b"test"

    def test_get_bins(self):
        test = get_bins(500, 100)
        assert len(test) == 5
        assert test[0].start == 0
        assert test[0].end == 100
        assert test[1].start == 100
        assert test[1].end == 200
        assert test[4].start == 400
        assert test[4].end == 500


class TestBedline(object):

    def test_str(self, test_bedline):
        assert str(test_bedline) == "chr1\t1\t100\t20"

    def test_bytes(self, test_bedline):
        assert bytes(test_bedline) == b"chr1\t1\t100\t20"

    def test_eq(self, test_bedline):
        assert test_bedline == test_bedline

    def test_neq(self, test_bedline):
        other_line = BedLine(b"chr2", 1, 100, 20)
        assert test_bedline != other_line

    def test_fromline(self, test_bedline):
        n = BedLine.fromline("chr1\t1\t100\t20")
        assert n == test_bedline
        n2 = BedLine.fromline(b"chr1\t1\t100\t20")
        assert n2 == test_bedline
        three_fields = BedLine.fromline("chr1\t1\t100")
        assert three_fields == BedLine(b"chr1", 1, 100, b"NA")

