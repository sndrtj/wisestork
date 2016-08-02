from collections import namedtuple

Bin = namedtuple("Bin", ["start", "end"])


class BedLine(namedtuple("BedLine", ["chromosome", "start", "end", "value"])):
    __slots__ = ()

    def __str__(self):
        return "{0}\t{1}\t{2}\t{3}".format(self.chromosome, self.start, self.end, self.value)


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