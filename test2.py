import argparse
import gzip

import numpy as np
import pysam

from utils import BedLine


def get_refbins_from_db(tbx_handle, bed_line):
    """
    Get reference bins from database
    :param tbx_handle: tabix_handle
    :param bed_line: BedLine namedtuple
    :return: list of BedLine namedtuples (may be empty)
    """
    record = [x for x in tbx_handle.fetch(str(bed_line.chromosome), bed_line.start, bed_line.end)][0]
    rows = [x.replace(",", "\t") for x in record.strip().split("\t")[3].split("|")]
    if rows == [""]:
        return []
    if rows == ['nan']:
        return []
    bedlines = [BedLine.fromline(x) for x in rows]
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-I", "--input", type=str, required=True, help="Input bgzipped and tabixxed bedgraph file")
    parser.add_argument("-D", "--database", type=str, required=True, help="Reference database")
    parser.add_argument("-O", "--output", type=str, required=True, help="Output file")

    args = parser.parse_args()
    db_handle = pysam.Tabixfile(args.database)

    ihandle = gzip.open(args.input)
    tbx = pysam.Tabixfile(args.input)
    ohandle = open(args.output, "w")
    for i, r in enumerate(ihandle):
        print i
        tmp = BedLine.fromline(r)
        lines = get_refbins_from_db(db_handle, tmp)
        reference_bins = []
        for line in lines:
            reference_bins += [BedLine.fromline(x) for x in tbx.fetch(str(line.chromosome), line.start, line.end)]
        z = get_z_score(tmp, reference_bins)
        n = BedLine(tmp.chromosome, tmp.start, tmp.end, z)
        ohandle.write(str(n) + "\n")

    ohandle.close()
    ihandle.close()
    tbx.close()
