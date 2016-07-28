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
    return [(s, e) for s, e in zip(starts, ends)]