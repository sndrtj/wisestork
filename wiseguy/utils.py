"""
wiseguy.utils
~~~~~~~~~~~~~

:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
:license: GPLv3
"""


from collections import namedtuple

Bin = namedtuple("Bin", ["start", "end"])


class BedLine(namedtuple("BedLine", ["chromosome", "start", "end", "value"])):
    __slots__ = ()

    def __str__(self):
        return "{0}\t{1}\t{2}\t{3}".format(self.chromosome, self.start, self.end, self.value)

    def __eq__(self, other):
        return self.chromosome == other.chromosome and self.start == other.start and self.end == other.end

    def __ne__(self, other):
        return not self.__eq__(other)

    @classmethod
    def fromline(cls, str):
        contents = str.split("\t")
        return cls(contents[0], int(contents[1]), int(contents[2]), float(contents[3]))


def attempt_integer(value):
    """
    Attempt to make an integer from a string
    return string when not possible
    :param value: the input string
    :return: the possibly-converted item
    """
    try:
        return int(value)
    except ValueError:
        return value


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
    return [Bin(s, e) for s, e in zip(starts, ends)]