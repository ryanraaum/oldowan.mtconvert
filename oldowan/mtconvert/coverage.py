from seq2sites import find_match_positions
from types import TupleType, IntType, ListType

# rCRSplus is an expanded rCRS sequence, which starts at position 15500,
#          then runs through the whole genome, then has the opening
#          1000 bases attached at the end again. This is used so that
#          query sequences that are not cut precisely at the canonical 
#          origin will still be analyzable.
# rCRSplus_positions maps the string indices of rCRSplus to their positions
#          in the reference sequence.
from oldowan.mtdna import rCRSplus, rCRSplus_positions

import itertools

WORD_SIZE = 15
FORCE_SPLIT = 50

def std_cmp(x,y):
    if x < y:
        return -1
    return 1

def hvr_cmp_f(x,y):
    if (x < 600 and y < 600) or (x > 15500 and y > 15500):
        return std_cmp(x,y)
    if x > 15500 and y < 600:
        return -1
    else:
        return 1

class Coverage(object):
    
    def __init__(self, *segments):
        self._sites = self._segments_to_sites(segments)
        self._str = self._make_str()

    def _segments_to_sites(segments):
        sites = []
        for s in segments:
            if type(s) in (TupleType,ListType):
                start = s[0]
                stop  = s[-1]
                if stop == 0:
                    stop = 16569
                if start > stop:
                    sites = sites + range(start,16570) + range(1,stop+1)
                else:
                    sites = sites + range(start,stop+1)
            if type(s) is IntType:
                sites.append(s)

        sites = list(set(sites))
        # little funkiness here to keep the order hvr1 -> hvr2 
        if all(x < 600 or x > 15500 for x in sites):
            sites.sort(cmp=hvr_cmp_f)
        else:
            sites.sort()
        return sites
    _segments_to_sites = staticmethod(_segments_to_sites)

    def __get_sites(self):
        return self._sites

    sites = property(fget=__get_sites)

    def __repr__(self):
        return self._str

    def __cmp__(self, other):
        if self._sites == other.sites:
            return 0
        elif self._sites < other.sites:
            return -1
        return 1

    def _make_str(self):
        s = []
        i = 0
        start = None
        stop  = None
        while i < len(self._sites):
            if start is None:
                start = self._sites[i]
            i += 1
            if i >= len(self._sites) or self._sites[i-1] + 1 != self._sites[i]:
                stop = self._sites[i-1]
                if start == stop:
                    s.append("%d" % start)
                else:
                    s.append("%d:%d" % (start,stop))
                start = None
                stop  = None
        return ','.join(s)

    def intersection(self, other):
        same_sites = list(set(self.sites).intersection(other.sites))
        return Coverage(*same_sites)

def calc_num_terminal_mismatches(matches):
    """(Internal) Count the number of -1 entries at the end of a list of numbers.

    These -1 correspond to mismatches in the sequence to sequence search.
    """
    if matches[-1] != -1:
        return 0
    neg_len = -1
    while matches[neg_len] == -1:
        neg_len -= 1
    return abs(neg_len+1)


def calc_num_opening_mismatches(matches):
    """(Internal) Count the number of -1 entries at the beginning of a list of numbers.

    These -1 correspond to mismatches in the sequence to sequence search.
    """
    if matches[0] != -1:
        return 0
    num = 1
    while matches[num] == -1:
        num += 1
    return num


def chunk(chunks, query, word_size=WORD_SIZE, force_split_at=FORCE_SPLIT):
    """(Internal) Calculate match positions of query against reference.

    Attempts to identify pasted together query sequences.
    """
    # first, find match positions
    matches = find_match_positions(query, rCRSplus, word_size)

    # if there are opening mismatches, need to fix the opening number
    num_opening_mismatches = calc_num_opening_mismatches(matches)
    if num_opening_mismatches > 0:
        matches[0] = matches[num_opening_mismatches] - num_opening_mismatches

    num_terminal_mismatches = calc_num_terminal_mismatches(matches)
    # find the end of the matching range
    # (could be the whole length if there are no terminal mismatches)
    if num_terminal_mismatches >= force_split_at:
        end_of_range = len(matches) - num_terminal_mismatches
        current_chunk = matches[:end_of_range]
        current_chunk.append(current_chunk[-1]+word_size-1)
        chunks.append(current_chunk)
        return chunk(chunks, query[-num_terminal_mismatches:], word_size)
    else:
        current_chunk = matches
        current_chunk.append(current_chunk[-(1+num_terminal_mismatches)]+num_terminal_mismatches+word_size-1)
        chunks.append(current_chunk)
        return chunks


def calc_coverage(query, word_size=WORD_SIZE, force_split_at=FORCE_SPLIT):
    chunks = []
    chunk(chunks, query, word_size, force_split_at)
    # chunks are numbered relative to rCRS_plus sequence
    # need to change that to canonical rCRS position numbering
    tr = rCRSplus_positions
    segments = list((tr[x[0]+1], tr[x[-1]+1]) for x in chunks)
    return Coverage(*segments)

