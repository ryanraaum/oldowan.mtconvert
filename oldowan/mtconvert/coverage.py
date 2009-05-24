from seq2sites import find_match_positions

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

class Segment(object):
    
    def __init__(self, start, stop):
        self.start = start
        self.stop = stop
        if stop == 0:
            # there is no zero.
            self.stop = 16569

    def __repr__(self):
        if self.start == self.stop:
            return "%d" % self.start
        return "%d:%d" % (self.start, self.stop)

    def to_site_list(self):
        if self.start > self.stop:
            return range(self.start, 16570) + range(1,self.stop+1)
        return range(self.start, self.stop+1)


def create_segment(item):
    """(Internal) Create a Segment instance from a tuple (start, stop).

    If argument supplied is already a segment instance, will pass it through
    unchanged.
    """
    if isinstance(item, Segment):
        return item
    return Segment(item[0], item[1])


class Coverage(object):
    
    def __init__(self, *segments):
        self.segments = list(create_segment(x) for x in segments)
        self._str = self._make_str()

    def __repr__(self):
        return self._str

    def __cmp__(self, other):
        my_sites = list(set(self.to_site_list()))
        my_sites.sort()
        other_sites = list(set(other.to_site_list()))
        other_sites.sort()
        if my_sites == other_sites:
            return 0
        elif my_sites < other_sites:
            return -1
        return 1

    def _make_str(self):
        return ','.join(list(str(s) for s in self.segments))

    def to_site_list(self):
        return list(itertools.chain(*list(x.to_site_list() for x in self.segments)))


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
    segments = list(Segment(tr[x[0]+1], tr[x[-1]+1]) for x in chunks)
    return Coverage(*segments)

