import re

from oldowan.polymorphism import Polymorphism

# rCRSlist is the rCRS sequence exploded into a list with biological numbering
from oldowan.mtdna import rCRSlist

RE_SITE = re.compile(r"([0-9]+)(.*)")

TRANSITIONS = { 'G':'C',
                'C':'G',
                'A':'T',
                'T':'A' }

def str2sites(a_string, add16k=False):
    """Extract a list of Polymorphisms from a string.

    """
    if isinstance(a_string, Polymorphism):
        return [a_string]

    # split the string by whitespace after replacing all commas and semi-colons
    # with whitespace
    entries = re.sub(r'[,;]', ' ', a_string).split()

    # if given 'rCRS', 'CRS', or the empty string '', there are no sites
    if (a_string.strip().upper() == 'RCRS' or 
        a_string.strip().upper() == 'CRS' or
        a_string.strip() == ''):
        return []

    # each position is found by a regex that looks for a number followed by
    # anything
    matches = list(RE_SITE.match(x) for x in entries)
    invalid_site_indices = list(pos for pos, val in enumerate(matches) if val == None)
    if len(invalid_site_indices) > 0:
        invalid_sites = list(entries[x] for x in invalid_site_indices)
        raise Exception('These sites are invalid: %s' % invalid_sites)
    positions = [RE_SITE.match(x).groups() for x in entries]

    # Populate the site list that we will return, making sure to fix the
    # edge cases of:
    # 1. no value for the position (we assume transition)
    # 2. a value larger than a single character, which could be a:
    #   a. an insertion, which could look like:
    #       '16129insA'
    #       '16129insAT'
    #       '16129+A'
    #       '16129+AT'
    #       '16129.1A'
    #       '16129.2T'
    #   b. a deletion, which could look like:
    #       '16129del'
    #       '16129delG'
    sites = []
    for k,v in positions:
        insert = 0
        k = int(k)
        if add16k:
            k = k+16000
        if v == '':
            v = TRANSITIONS[rCRSlist[k]]
        elif v == 'd':
            v = '-'
        elif len(v) > 1:
            if v.upper().startswith('DEL'):
                v = '-'
            elif v.startswith('.'):
                # this is probably an insertion. i.e. 16129.1A
                match = re.match(r'\.([1-9][0-9]*)([ACGTURYMKSWBDHVN]+)',v)
                if match:
                    insert, v = match.groups()
                    insert = int(insert)
            elif v.upper().startswith('INS'):
                # this is probably an insertion. i.e. 16129insA
                v = v.upper()
                match = re.match(r'INS([ACGTURYMKSWBDHVN]+)',v)
                if match:
                    (v,) = match.groups()
                    insert = 1
            elif v.startswith('+'):
                # this is probably an insertion. i.e. 16129+A
                insert = True
                match = re.match(r'\+([ACGTURYMKSWBDHVN]+)',v)
                if match:
                    (v,) = match.groups()
                    insert = 1

        if v == '-' or insert == 0:
            sites.append(Polymorphism(k,0,v))
        elif insert > 0:
            if len(v) > 1:
                for i in range(len(v)):
                    sites.append(Polymorphism(k,i+1,v[i]))
            else:
                sites.append(Polymorphism(k,insert,v))

    return sites

