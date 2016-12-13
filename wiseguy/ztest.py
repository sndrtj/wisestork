"""
wiseguy.ztest
~~~~~~~~~~~~~

:copyright: (c) 2013 Roy Straver
:copyright: (c) VU University Medical Center
:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
:license: GPLv3
"""

import argparse
import gzip
import sys

import numpy as np
import pysam

from .utils import BedLine


def get_refbins_from_db(tbx_handle, bed_line):
    """
    Get reference bins from database
    :param tbx_handle: tabix_handle
    :param bed_line: BedLine namedtuple
    :return: list of BedLine namedtuples (may be empty)
    """
    record = [x for x in tbx_handle.fetch(bed_line.chromosome.decode(), bed_line.start, bed_line.end)][0]
    rows = [x.replace(",", "\t") for x in record.strip().split("\t")[3].split("|")]
    if rows == [""]:
        return []
    if rows == ['nan']:
        return []
    bedlines = [BedLine.fromline(x.encode()) for x in rows]
    return bedlines


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


def ztest(input_path, output_path, database_path, counter_interval=1000):
    """
    Calculate z scores from bed file and database bed file
    :param input_path: query bed file path
    :param output_path: output bed file path
    :param database_path: database file path
    :param counter_interval: Interval for counter updates (default: 1000)
    :return: -
    """
    db_handle = pysam.Tabixfile(database_path)
    ihandle = gzip.open(input_path)
    tbx = pysam.Tabixfile(input_path)
    ohandle = open(output_path, "wb")
    for i, r in enumerate(ihandle):
        if i % counter_interval == 0:
            print("Processed {i} records".format(i=i), file=sys.stderr)
        tmp = BedLine.fromline(r)
        lines = get_refbins_from_db(db_handle, tmp)
        reference_bins = []
        for line in lines:
            reference_bins += [BedLine.fromline(x.encode()) for x in tbx.fetch(line.chromosome.decode(), line.start, line.end)]
        z = get_z_score(tmp, reference_bins)
        n = BedLine(tmp.chromosome.decode(), tmp.start, tmp.end, z)
        ohandle.write(bytes(n) + b"\n")

    ohandle.close()
    ihandle.close()
    tbx.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-I", "--input", type=str, required=True, help="Input bgzipped and tabixxed bedgraph file")
    parser.add_argument("-D", "--database", type=str, required=True, help="Reference database")
    parser.add_argument("-O", "--output", type=str, required=True, help="Output file")

    args = parser.parse_args()