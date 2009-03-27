import re
from types import StringType, ListType, DictType, IntType

from oldowan.polymorphism import Polymorphism
from str2sites import str2sites

# rCRSlist is the rCRS sequence exploded into a list with biological numbering
from oldowan.mtdna import rCRSlist

from oldowan.mtdna import HVR1_indices
from oldowan.mtdna import HVR2_indices
from oldowan.mtdna import HVR1and2_indices
from oldowan.mtdna import HVR1to2_indices
from oldowan.mtdna import coding_indices
from oldowan.mtdna import all_indices

REGIONS = { 'HVR1'    : HVR1_indices, 
            'HVR2'    : HVR2_indices,
            'HVR1AND2': HVR1and2_indices,
            'HVR1TO2' : HVR1to2_indices,
            'CODING'  : coding_indices,
            'ALL'     : all_indices }

VALID_REGIONS = REGIONS.keys()

def flatten(item):
    if type(item) == ListType:
        if len(item) > 1:
            return flatten(item[0]) + flatten(item[1:])
        elif len(item) == 1:
            return flatten(item[0])
        else:
            return []
    else:
        return [item]


def sites2seq(sites, region='HVR1', add16k=False):
    """
    Translate sites to sequence.

    The 'sites' variable can be provided as: 
    1. a string of positions and their values:
        '16129A 16189C 16223T'
    2. a list of strings:
        ['16129A', '16189C', '16223T']
    2. a list of Polymorphisms:
        [Polymorphism(16129,0,'A'), Polymorphism(16189,0,'C'), ... ]

    For the revised Cambridge Reference Sequence, either submit the empty 
    string ('') or 'rCRS' as the sites variable.

    The options for the 'region' argument are:
    1. 'HVR1'
    2. 'HVR2'
    3. 'HVR1and2'
    4. 'HVR1to2'
    5. 'coding'
    6. 'all'
    7. an array of sites. These site numbers are assumed to start at '1'
       following the usual mtDNA site numbering convention. (Not the starts
       at '0' computer science convention.)

    I am following the EMPOP Mitochondrial DNA Control Region Database 
    definition of HVR1 (HVS-I) and HVR2 (HVS-II).

    HVR1 consists of bases [16024-16365] of the revised Cambridge Reference
    Sequence (rCRS), which is available in GenBank under accession number
    AC_000021.2 (gi: 115315570)

    HVR2 consists of bases [73-340] of the rCRS.

    HVR1and2 is HVR1 and HVR2 concatenated.

    HVR1to2 runs from the start of HVR1 around to the end of HVR2.


    Calling sites2seq with an empty string or array will return 'region' from
    the rCRS. The default 'region' is HVR1.

    """
    region_type = 'invalid'
    if type(region) == ListType:
        region_type = 'list'
    elif callable(getattr(region, 'upper')) and region.upper() in VALID_REGIONS:
        region_type = 'string'
        region = region.upper()
    if region_type == 'invalid':
        raise Exception('"region" argument "%s" is invalid' % region)

    # the rCRSlist sequence has had an opening spacer added 
    # so that indexing starts at 1
    working = list(rCRSlist) 

    if (type(sites) == StringType):
        if sites.upper() == 'RCRS':
            sites = ''
        sites = str2sites(sites,add16k=add16k)
    elif (type(sites) == ListType):
        sites = flatten(list(str2sites(x,add16k=add16k) for x in sites))
    else:
        raise Exception('"sites" argument "%s" is not acceptable' % type(sites))

    # need to deal with insertions separately so that the index numbering
    # can be maintained
    insertions = []
    
    for poly in sites:
        if poly.position < 1 or poly.position > 16569:
            msg = "'%s' is not in the human mtDNA sequence range" % poly
            raise Exception, msg
        elif poly.insert > 0:
            insertions.append(poly)
        else:
            working[poly.position] = poly.value

    requested_positions = None
    if region_type == 'string':
        requested_positions = REGIONS[region]
    else:
        requested_positions = region

    # as insertions are added to the working sequence, we have to modify the
    # requested_positions to account for changes in the length of the sequence
    if len(insertions) > 0:
        # working from highest position to lowest means that we don't have to
        # keep updating the position of the insertion
        insertions.sort(reverse=True)
        for current in insertions:
            increment = lambda x: x + 1 if x > current.position else x
            requested_positions = list(increment(x) for x in requested_positions)
            # the list completion properly increments all positions higher than
            # the insertion position, but ends up skipping the point of
            # insertion altogether, so add that back in to the list
            insert_index = requested_positions.index(current.position)+1
            requested_positions.insert(insert_index, current.position+1)
            working.insert(current.position+1, current.value)
        
    result = list(working[x] for x in requested_positions)
    return ''.join(result)

