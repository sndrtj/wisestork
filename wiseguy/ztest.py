"""
wiseguy.ztest
~~~~~~~~~~~~~

:copyright: (c) 2013 Roy Straver
:copyright: (c) VU University Medical Center
:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
:license: GPLv3
"""

import gzip
from math import isnan
import sys

import numpy as np
import progressbar

from .utils import BedLine, utf8


def get_z_score(bin, reference_bins):
    """
    Get z score of for this bin
    :param bin: bin
    :param reference_bins: reference bins
    :return: float
    """
    if len(reference_bins) == 0:
        return np.nan
    average = np.average([x.value for x in reference_bins])
    stddev = np.std([x.value for x in reference_bins])

    if stddev == 0:
        return 0
    Z = (bin.value - average) / stddev
    return Z


def ztest(input_path, output_path, database_path):
    """
    Calculate z scores from gzipped bed file and database bed file
    :param input_path: query bed file path
    :param output_path: output bed file path
    :param database_path: database file path
    :return: -
    """
    ohandle = open(output_path, "wb")
    refdict = build_reference_index(database_path)
    with gzip.open(input_path) as ihandle:
        bedlines = [BedLine.fromline(x) for x in ihandle]

    if not len(bedlines) == len(refdict):
        raise ValueError("Reference and query bed files "
                         "are of different size!")

    print("Calculating Z-scores")
    bar = progressbar.ProgressBar(max_value=len(bedlines))

    for i, x in enumerate(bedlines):
        bar.update(i)
        key = create_key(x)
        reference_bins = [bedlines[i] for i in refdict[key]]
        z = get_z_score(x, reference_bins)
        n = BedLine(x.chromosome.decode(), x.start, x.end, z)
        ohandle.write(bytes(n) + b"\n")
    bar.finish()
    ohandle.close()


def create_key(bedline):
    return b"|".join(map(utf8, [bedline.chromosome, bedline.start, bedline.end]))


def build_reference_index(database_path):
    """
    Build a reference index from a database path.

    The index is a dictionary of type :
    {key: [list of indices]}

    where the key is a bed region, and the indices are indices of similar
    bins. Naturally, this means that query bed files _must_ be of
    _exactly_ the same size as the database, and the sorting order _must_
    also be identical.

    :param database_path:
    :return: index dictionary
    """
    print("Building index", file=sys.stderr)
    n = 0
    d = {}
    with gzip.open(database_path) as db_handle:
        for i, line in enumerate(db_handle):
            b = BedLine.fromline(line)
            key = create_key(b)
            d[key] = i
            n = i + 1
    refdict = {}
    with gzip.open(database_path) as db_handle2:
        bar = progressbar.ProgressBar(max_value=n)
        for o, l in enumerate(db_handle2):
            bar.update(o)
            b = BedLine.fromline(l)
            key2 = create_key(b)
            if isinstance(b.value, float) and isnan(b.value):
                refdict[key2] = []
            else:
                it = [BedLine.fromline(x, b",") for x in b.value.split(b"|")]
                keys = list(map(create_key, it))
                idxs = [d[x] for x in keys]
                refdict[key2] = idxs
        bar.finish()

    return refdict