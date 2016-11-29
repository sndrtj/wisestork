"""
wiseguy.wiseguy
~~~~~~~~~~~~~~~~~~~~~

:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
"""

import click

from functools import wraps

from . import __version__
from .count import count
from .gc_correct import gc_correct

shared_options = [
    click.option("--binsize", "-B", type=click.IntRange(0, None), default=50000,
                 help="Bin size to use. Default = 50000"),
    click.option("--reference", "-R", type=click.Path(exists=True), required=True,
                 help="Path to reference fasta"),
    click.version_option(version=__version__)
]


def generic_option(options):
    """
    Decorator to add generic options to Click CLI's
    The group parent should NOT be decorated with this decorator
    :param options: list of click.option
    :return: decorated function
    """
    def __generic_option(func):
        for option in reversed(options):
            func = option(func)
        return func
    return __generic_option


@click.command(short_help="Count coverages")
@generic_option(shared_options)
@click.option("--output", "-O", type=click.Path(), required=True, help="Path to output BED file")
@click.option("--input", "-I", type=click.Path(exists=True), required=True, help="Path to input BAM file")
def count_cli(**kwargs):
    """
    Take a BAM file, and calculate the number of reads per bin.
    It will output a BED file (with 0-based positions) of regions and
    associated number of reads per bin.

    \b
    Your BAM file _must_ be indexed and _must_ contain chromosome
    lengths and names in the header.
    """
    input = kwargs.get("input", None)
    output = kwargs.get("output", None)
    binsize = kwargs.get("binsize", 50000)
    reference = kwargs.get("reference", None)
    count(input=input, output=output, binsize=binsize, reference=reference)


@click.command(short_help="GC correct")
@generic_option(shared_options)
@click.option("--output", "-O", type=click.Path(), required=True, help="Path to output BED file")
@click.option("--input", "-I", type=click.Path(exists=True), required=True, help="Path to input BED file")
@click.option("--frac-n", "-n", type=click.FLOAT, default=0.1, help="Maximum fraction of N-bases per bin. "
                                                                    "Default = 0.1")
@click.option("--frac-r", "-r", type=click.FLOAT, default=0.0001, help="Minimum fraction of reads per bin. "
                                                                       "Default = 0.0001")
@click.option("--iter", "-t", type=click.INT, default=3, help="Number of iterations for LOWESS function. Default = 3")
@click.option("--frac-lowess", "-l", type=click.FLOAT, default=0.1, help="Fraction of data to use for LOWESS function. "
                                                                         "Default = 0.1")
def gcc_cli(**kwargs):
    """
    GC-correct a BED file containing counts per region
    """
    input_path = kwargs.get("input", None)
    output = kwargs.get("output", None)
    reference = kwargs.get("reference", None)
    frac_n = kwargs.get("frac_n", 0.1)
    frac_r = kwargs.get("frac_r", 0.0001)
    iter = kwargs.get("iter", 3)
    frac_lowess = kwargs.get("frac_lowess", 0.1)
    gc_correct(input=input_path, output=output, reference=reference,
               frac_r=frac_r, frac_n=frac_n, iter=iter, frac_lowess=frac_lowess)


@click.command(short_help="Calculte Z-scores")
@generic_option(shared_options)
def zscore_cli(**kwargs):
    pass

@click.command(short_help="Create new reference")
@generic_option(shared_options)
def newref_cli(**kwargs):
    pass


@click.group()
@click.version_option()
def cli(**kwargs):
    """
    Discover CNVs from BAM files.

    \b
    A typical workflow first extracts regions from a BAM file
    The resulting BED tracks must then be GC-corrected.
    Using a reference track of region similarity,
    One can then calculate Z-scores for every region.

    \b
    The following sub-commands are supported:
     - count: count coverage per bin
     - gc-correct: GC-correct bins
     - zscore: calculate Z-scores
     - newref: Generate a new reference dictionary of bin similarities

    """
    pass

def main():
    cli.add_command(count_cli, "count")
    cli.add_command(gcc_cli, "gc-correct")
    cli.add_command(zscore_cli, "zscore")
    cli.add_command(newref_cli, "newref")
    cli()

if __name__ == '__main__':
    main()
