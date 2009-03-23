import re

from types import StringType,ListType

import biopython.pairwise2

from oldowan.polymorphism import Polymorphism

from reduction_funcs import prefer_known_substitutions
from reduction_funcs import prefer_fewer
from reduction_funcs import prefer_indels_at_end
from reduction_funcs import prefer_multi_inserts
from reduction_funcs import prefer_insertions_at_309_and_315
from reduction_funcs import prefer_315_insertion_over_double_310_insertion
from reduction_funcs import prefer_indels_over_ts_over_tv 
from reduction_funcs import prefer_end_of_sorted_list
from reduction_funcs import prefer_309del_315ins_over_309T_310C
from reduction_funcs import prefer_95ins_97del_over_96T_97C

# rCRSplus is an expanded rCRS sequence, which starts at position 15500,
#          then runs through the whole genome, then has the opening
#          1000 bases attached at the end again. This is used so that
#          query sequences that are not cut precisely at the canonical 
#          origin will still be analyzable.
# rCRSplus_positions maps the string indices of rCRSplus to their positions
#          in the reference sequence.
from oldowan.mtdna import rCRSplus, rCRSplus_positions

class MatchingRange(object):
    
    def __init__(self, query_start, target_start):
        self.query_start = query_start
        self.target_start = target_start
        self.query_end = None
        self.target_end = None
        
    def __str__(self):
        return '[(%s, %s),(%s,%s)]' % (self.query_start, 
                                       self.query_end, 
                                       self.target_start, 
                                       self.target_end)
    
    def __repr__(self):
        return 'MatchingRange: %s' % str(self)

    def query_slice(self):
        return slice(self.query_start, self.query_end)

    def target_slice(self):
        return slice(self.target_start, self.target_end)

    def intersect(self, other):
        other_query_range  = range(other.query_start, other.query_end)
        other_target_range = range(other.target_start, other.target_end)
        if ((self.query_start  in other_query_range  or
             self.query_end    in other_query_range )  and
            (self.target_start in other_target_range or
             self.target_end   in other_target_range)  ):
            return True
        return False

    def merge_with(self, other):
        self.query_start  = min([self.query_start, other.query_start])
        self.query_end    = max([self.query_end, other.query_end])
        self.target_start = min([self.target_start, other.target_start])
        self.target_end   = max([self.target_end, other.target_end])

def find_match_positions(query, reference, word_size=15):

    WS = word_size
    s1 = query
    s2 = reference
    
    word_starts = range(0, len(s1)-WS+1) 
    matches = []
    start = None
    for i in word_starts:
        if start is None:
            pos = s2.find(s1[i:i+WS])
        else:
            pos = s2.find(s1[i:i+WS], start, start+100)
        if pos != -1:
            start = pos
        matches.append(pos)

    return matches


def align(query, reference, word_size=15, mismatch_cutoff=0.7):
    """
    Align two similar sequences.
    
    """
    WS = word_size
    s1 = query
    s2 = reference
    
    matches = find_match_positions(s1, s2, WS)

    if -1 not in matches:
        return []
    
    if matches.count(-1)/float(len(matches)) > mismatch_cutoff:
        raise Exception('sequences do not match')

    mismatches = []
    mismatch = None
    for position, value in enumerate(matches):
        if value == -1:
            if mismatch is None:
                if position == 0:
                    # the query starts at 0, but the target doesn't have to
                    # so, figure out where the target start is
                    query_start = 0
                    count = 0
                    while matches[count] == -1:
                        count += 1
                    target_start = matches[count] - count 
                    mismatch = MatchingRange(query_start, target_start)
                else:
                    mismatch = MatchingRange(position, matches[position-1]+1)
        elif mismatch is not None:
            mismatch.query_end = position+WS
            mismatch.target_end = value+WS
            mismatches.append(mismatch)
            mismatch = None

    # if there is a mismatch at the end of the query, we will
    # never have set the proper end of the query and target
    if mismatch is not None:
        mismatch.query_end = len(query)
        mismatch.target_end = (mismatch.target_start + 
                               mismatch.query_end - 
                               mismatch.query_start)
        mismatches.append(mismatch)
        mismatch = None

    polymorphisms = []

    merged_matches = []
    for mm in mismatches:
        merged = False
        for nm in merged_matches:
            if mm.intersect(nm):
                nm.merge_with(mm)
                merged = True
        if not merged:
            merged_matches.append(mm)
    
    for mm in merged_matches:
        slice1 = s1[mm.query_slice()]
        slice2 = s2[mm.target_slice()]
        # align2 args (seq1, seq2, match_score, mismatch_penalty, gap_penalty, gap_extension_penalty)
        alignments = biopython.pairwise2.align.globalms(slice1, slice2, 
                                                        3, -1, -3, -1)
        aln_polymorphisms = []        
        for alignment in alignments:
            this_aln_polymorphisms = []
            qry_aln, ref_aln = alignment[0], alignment[1]
            abs_pos = mm.target_start - 1
            insert = 0
            for pos,val in enumerate(ref_aln):
                if val != qry_aln[pos]:
                    if val == '-':
                        insert += 1
                    else:
                        abs_pos += 1
                        insert = 0
                    new_poly = Polymorphism(abs_pos,insert,qry_aln[pos]) 
                    this_aln_polymorphisms.append(new_poly)
                else:
                    abs_pos += 1
                    insert = 0
            aln_polymorphisms.append(this_aln_polymorphisms)
        
        polymorphisms.append(aln_polymorphisms)

    return polymorphisms


def seq2sites(seq, word_size=15, mismatch_cutoff=0.7, ambig_cutoff=10):
    
    if seq.count('N') > ambig_cutoff:
        raise Exception, "too many N's in submitted sequence"

    # remove whitespace
    seq = re.sub(r'\s', '', seq.upper())

    polymorphisms = align(seq, rCRSplus, word_size, mismatch_cutoff)

    for mismatch_block in polymorphisms:
        for alternate_alignment in mismatch_block:
            for pos,site in enumerate(alternate_alignment):
                # +1 to move from zero-based to one-based counting
                site.position = rCRSplus_positions[site.position] + 1

    reduction_funcs = [ prefer_known_substitutions,
                        prefer_insertions_at_309_and_315,
                        prefer_315_insertion_over_double_310_insertion,
                        prefer_309del_315ins_over_309T_310C,
                        prefer_95ins_97del_over_96T_97C,
                        prefer_fewer,
                        prefer_multi_inserts,
                        prefer_indels_at_end,
                        prefer_indels_over_ts_over_tv,
                        prefer_end_of_sorted_list ]

    polys = []
    for block in polymorphisms:
        for f in reduction_funcs:
            block = f(block)
        polys.append(block)

    # elminate a lot of the excess nesting
    unnested = []
    for mismatch_block in polys:
        if len(mismatch_block) == 1:
            for site in mismatch_block[0]:
                unnested.append(site)
        else:
            unnested.append(mismatch_block)
    return unnested

