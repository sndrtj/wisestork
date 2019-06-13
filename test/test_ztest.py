from collections import namedtuple
import random

import pytest

from wisestork.ztest import *
from wisestork.utils import BedLine

ValueObject = namedtuple("ValueObject", ['value'])  # little helper object


class TestFunctions:

    def test_create_key(self):
        bline = BedLine(b"chr1", 1, 100, 0)
        bline2 = BedLine(b"chr1", b"1", b"100", 0)
        assert create_key(bline) == b"chr1|1|100"
        assert create_key(bline2) == b"chr1|1|100"

    def test_build_reference_index(self):
        pass

    def test_get_z_score(self):
        assert isnan(get_z_score(None, []))
        objects = [ValueObject(random.normalvariate(100, 20)) for _ in range(2000)]
        zero_test = ValueObject(100)
        assert -0.5 < get_z_score(zero_test, objects) < 0.5
        minus_test = ValueObject(20)
        assert -4.5 < get_z_score(minus_test, objects) < -3.5
        positive_test = ValueObject(160)
        assert 2.5 < get_z_score(positive_test, objects) < 3.5
        zero_std = [ValueObject(random.normalvariate(100, 0)) for _ in range(2000)]
        assert isnan(get_z_score(zero_test, zero_std))

    def test_ztest(self):
        pass
