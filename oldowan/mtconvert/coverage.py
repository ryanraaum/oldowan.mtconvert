from seq2sites import find_match_positions

# rCRSplus is an expanded rCRS sequence, which starts at position 15500,
#          then runs through the whole genome, then has the opening
#          1000 bases attached at the end again. This is used so that
#          query sequences that are not cut precisely at the canonical 
#          origin will still be analyzable.
# rCRSplus_positions maps the string indices of rCRSplus to their positions
#          in the reference sequence.
from oldowan.mtdna import rCRSplus, rCRSplus_positions

WORD_SIZE = 15
FORCE_SPLIT = 50

def calc_num_terminal_mismatches(matches):
    if matches[-1] != -1:
        return 0
    neg_len = -1
    while matches[neg_len] == -1:
        neg_len -= 1
    return abs(neg_len+1)


def chunk(chunks, query, word_size=WORD_SIZE, force_split_at=FORCE_SPLIT):
    matches = find_match_positions(query, rCRSplus, word_size)
    num_terminal_mismatches = calc_num_terminal_mismatches(matches)
    end_of_range = len(matches) - num_terminal_mismatches
    chunks.append(matches[:end_of_range])
    if num_terminal_mismatches < force_split_at:
        return chunks
    else:
        return chunk(chunks, query[-num_terminal_mismatches:], word_size)


def find_range(query, word_size=WORD_SIZE, force_split_at=FORCE_SPLIT):
    chunks = []
    chunk(chunks, query, word_size, force_split_at)
    print chunks
    return chunks
