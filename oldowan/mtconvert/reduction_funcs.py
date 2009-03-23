from oldowan.polymorphism import Polymorphism
from mitomap_data import substitutions

# COMMON ABBREVIATIONS USED HERE
#
# mb: mismatch block, a list of lists
#     each sublist holds the Polymorphisms for an alternate alignment
#

#
# UTILITY FUNCTIONS
#

def cmp_alternate_alignments(a,b):
    """Compare function for lists of Polymorphisms.

    Comparison is based upon Polymorphism positions.
    """
    a.sort()
    b.sort()
    for i in range(len(a)-1,-1,-1):
        if a[i].position != b[i].position:
            return cmp(a[i].position, b[i].position)
        if a[i].insert != b[i].insert:
            return cmp(a[i].insert, b[i].insert)
    return 0

#
# REDUCTION FUNCTIONS
#

def prefer_fewer(mb):
    """Prefer fewer polymorphisms."""
    mb.sort(cmp=lambda x,y: cmp(len(x),len(y)))
    mb = list(x for x in mb if len(x) == len(mb[0]))
    return mb


def prefer_known_substitutions(mb):
    """Prefer substitutions previously seen to new substitutions.

    Previously seen is based upon the listing from www.mitomap.org.
    Each alternate alignment receives a score calculated based on
    how many of the substitutions in that alignment are in the
    mitomap variants list: 1 point for transversions and 2 points 
    for transitions.
    """
    if len(mb) > 1:
        scores = [0] * len(mb)
        # pos: position
        # aa: alternate alignment
        for pos, aa in enumerate(mb):
            for variant in aa:
                if variant.is_substitution():
                    if str(variant) in substitutions:
                        # larger scores are "better"
                        if variant.is_transition():
                            scores[pos] += 2
                        else:
                            scores[pos] += 1
        if max(scores) > 0:
            # lsi: indices of lower scoring alternate alignments
            lsi = list(x for x in range(len(mb)) if scores[x] < max(scores))
            # remove low scoring alignments from the mismatch block
            # in reverse order so as to not mess up the preceding indices
            lsi.sort(reverse=True)
            for i in lsi:
                mb.pop(i)

    return mb


def prefer_multi_inserts(mb):
    """Prefer several insertions at one site over separate insertions.

    When there is a choice between multiple inserts at one position
    or multiple inserts spread over several positions, choose
    multiple inserts at the same position.
    """
    # position_counts: number of insertion sites required in 
    #                  each alternate alignment
    position_counts = {}
    # pos: position
    # aa: alternate alignment
    for pos,aa in enumerate(mb):
        inserts = list(x for x in aa if x.insert > 0)
        positions = len(set(list(x.position for x in inserts)))
        position_counts[pos] = positions

    # if there are no insertions, this rule does not
    # apply, and there will be no position counts > 0
    if len(position_counts.values()) > 0:
        # minpos: minimum number of insertion sites
        minpos = min(position_counts.values())
    
        # only continue if there is any variation in number of
        # insertion positions required.
        if minpos != max(position_counts.values()):
            new_mb = []
            for i in range(len(mb)):
                if position_counts[i] == minpos:
                    new_mb.append(mb[i])
            mb = new_mb

    return mb


def prefer_indels_at_end(mb):
    """Prefer alignments with insertions and deletions moved 3'.

    For example, if an insertion could be placed in more than one position,
    choose the most 3' (3 prime) position. Here, the reference sequence has 
    two C's in a row, but the submitted sequence has three, so there are 3 
    possible alignments:

        Reference  : GA-CCTT
        Alignment 1: GACCCTT

        Reference  : GAC-CTT
        Alignment 2: GACCCTT

        Reference  : GACC-TT
        Alignment 3: GACCCTT

    We should choose Alignment 3 as this alternative has the 
    insertion in the most 3' possible position.
    """
    # don't waste time if there are no alternatives
    if len(mb) <= 1:
        return mb

    # indel_set_flag: set to True if the first alternate alignment
    #                 includes at least one indel. Falsified if later
    #                 alignments have a different number of indels than
    #                 the first. This rule will then not be applied.
    indel_set_flag = False

    # indel_set_len: the number of indels in the first alternative
    #                alignment. The others must have the same number.
    indel_set_len = None

    # indel_sets: the lists of indels in each alternate alignment is
    #             collected here
    indel_sets = []

    # aa: alternate alignment
    for aa in mb:

        indels = list(x for x in aa if x.insert > 0 or x.value == '-')

        # initialize the size of the insertion set
        if indel_set_len is None and len(indels) > 0:
            indel_set_flag = True
            indel_set_len = len(indels)

        # extract out sets of insertions
        if indel_set_flag and len(indels) == indel_set_len:
            indels.sort()
            indel_sets.append(list(x for x in indels))
        else:
            indel_set_flag = False

    # only apply the rule if we have found alternatives with indels
    if indel_set_len > 0 and indel_set_flag:
        indel_sets.sort(cmp=cmp_alternate_alignments)

        # full_sets: the whole alternative alignment, including substitutions
        full_sets = []
        for aa in mb:
            if indel_sets[-1][-1] in aa:
                full_sets.append(aa)

        full_sets.sort(cmp=cmp_alternate_alignments)

        return [full_sets[-1]]
        #return list(full_sets[-1] for i in range(len(mb)))
        
    return mb


def prefer_insertions_at_309_and_315(mb):
    """Prefer alternatives that include 309.1C or 315.1C over others.

    There are two multi-C runs at the beginning of the 300's and by
    convention, any insert in one of those runs is pushed to the end.
    Mostly, the other rules will pick these, but in some circumstances
    with other substitutions or deletions in this area, these won't get
    picked - although we want them to. Thus, this special case preference.
    """
    special_cases = ['309.1C', '315.1C']
    # mb: mismatch block 
    if len(mb) > 1:
        scores = [0] * len(mb)
        # pos: position
        # aa: alternate alignment
        for pos, aa in enumerate(mb):
            for variant in aa:
                if str(variant) in special_cases:
                    scores[pos] += 1
        if max(scores) > 0:
            # lsi: indices of lower scoring alternate alignments
            lsi = list(x for x in range(len(mb)) if scores[x] < max(scores))
            # remove low scoring alignments from the mismatch block
            # in reverse order so as to not mess up the preceding indices
            lsi.sort(reverse=True)
            for i in lsi:
                mb.pop(i)

    return mb


def prefer_315_insertion_over_double_310_insertion(mb):
    """Solely to (partially) fulfill Wilson et al. 2002 example 19.

    Wilson et al's solution to their example 19 fails to adhere to
    their rules. They suggest that the proper alignment should be:

        REF: CCT-CCCCC-GCT
        OBS: CCTTCCCCCCGCT

    However, their rule 3 is that insertions and deletions should
    be combined where the same number of differences to the reference
    sequence is maintained.

    Thus, following their rules, the alignment should be:

        REF: CCT--CCCCCGCT
        OBS: CCTTCCCCCCGCT

    However, Wilson et al 2002 have an unwritten rule 0:

        Rule 0: Prefer known variants.

    Thus, because the insertion placement rules will almost always
    put an insertion in the second run of C's above at 315.1C, as 
    in their preferred alignment, we will do the same here.

    However, any universal application of a rule that will do this
    breaks a lot of other desired behavior, thus, this special pleading.
    """
    before = ['310.1T', '310.2C']
    if len(mb) == 1 and len(mb[0]) == 2:
        var1, var2 = mb[0][0], mb[0][1]
        if str(var2) in before and str(var2) in before:
            var1 = Polymorphism(310,1,'T','-')
            var2 = Polymorphism(315,1,'C','-')
            return [[var1,var2]]
    return mb


def prefer_309del_315ins_over_309T_310C(mb):
    """Fulfill Wilson et al. 2002 example 18.

    Pairwise alignment parameters must favor substitutions over
    indels for most outcomes to be what is expected. However, there
    are situations, like this, where that mismatch between what must
    be done algorithmically conflicts with the selection rules (which
    state that indels should be preferred over substitutions). Thus,
    we have this special case.

    Without this special case, the solution to Wilson et al's example
    18 is 309T, 310C:

        REF: AAACCCCCCCTCCCCCGC
        OBS: AAACCCCCCTCCCCCCGC

    But, following the rule, that indels are preferred over substitutions,
    Wilson et al offer 309d, 315.1C:

        REF: AAACCCCCCCTCCCCC-GC
        OBS: AAACCCCCC-TCCCCCCGC

    However, their rule 3 is that insertions and deletions should
    be preferred where the same number of differences to the reference
    sequence is maintained.

    Therefore, this rule replaces 309T, 310C with 309d, 3151.C

    Any universal application of a rule that will do this
    breaks a lot of other desired behavior, thus, this special pleading.
    """
    before = ['309T', '310C']
    if len(mb) == 1 and len(mb[0]) == 2:
        var1, var2 = mb[0][0], mb[0][1]
        if str(var2) in before and str(var2) in before:
            var1 = Polymorphism(309,0,'-','C')
            var2 = Polymorphism(315,1,'C','-')
            return [[var1,var2]]
    return mb


def prefer_95ins_97del_over_96T_97C(mb):
    """Fulfill Wilson et al. 2002 example 22.

    Pairwise alignment parameters must favor substitutions over
    indels for most outcomes to be what is expected. However, there
    are situations, like this, where that mismatch between what must
    be done algorithmically conflicts with the selection rules (which
    state that indels should be preferred over substitutions). Thus,
    we have this special case.

    Without this special case, the solution to Wilson et al's example
    18 is 96T, 97C:

        REF: AGACGCTGGAGCC
        OBS: AGATCCTGGAGCC

    But, following the rule, that indels are preferred over substitutions,
    Wilson et al offer 95.1T, 97d:

        REF: AGA-CGCTGGAGCC
        OBS: AGATC-CTGGAGCC

    However, their rule 3 is that insertions and deletions should
    be preferred where the same number of differences to the reference
    sequence is maintained.

    Therefore, this rule replaces 96T, 97C with 95.1T,97d 

    Any universal application of a rule that will do this
    breaks a lot of other desired behavior, thus, this special pleading.
    """
    if len(mb) == 1:
        mb_str = list(str(x) for x in mb[0])
        if '96T' in mb_str and '97C' in mb_str:
            var1_ind = mb_str.index('96T')
            var2_ind = mb_str.index('97C')
            var1 = Polymorphism(95,1,'T','-')
            var2 = Polymorphism(97,0,'-','G')
            mb[0][var1_ind] = var1
            mb[0][var2_ind] = var2
    return mb


def prefer_16183C_16188T_16189C_over_16183d_16193insC(mb):
    """16183C and 16189C are common substitutions that are lost sometimes.

    This special case breaks the rules (more changes than necessary),
    but keeps common substitutions 16183C and 16189C, which is observed
    in 655 out of 5049 haplotypes in the 2007 Genographic haplotype listing.

    In this particular example, there is the 16188T substitution, which 
    leads the algorithm to choose 16183d, 16193.1C.
    """
    for aa in mb:
        aa_str = list(str(x) for x in aa)
        a = '16183d'
        b = '16193.1C'
        if a in aa_str and b in aa_str:
            var1_ind = aa_str.index(a)
            var2_ind = aa_str.index(b)
            indices = [var1_ind, var2_ind]
            indices.sort(reverse=True)
            for i in indices:
                aa.pop(i)
            aa.append(Polymorphism(16183,0,'C','A'))
            aa.append(Polymorphism(16188,0,'T','C'))
            aa.append(Polymorphism(16189,0,'C','T'))
            aa.sort() 
    return mb



def prefer_indels_over_ts_over_tv(mb):
    """Prefer indels over transitions over transversions.

    Following reccomendation 2 of Wilson et al. 2002:
       
       If there is more than one way to maintain the same number of 
       differences with respect to the reference sequence, differences
       should be prioritized in the following manner:
       (a) insertions/deletions (indels)
       (b) transitions
       (c) transversions

    """
    if len(mb) > 1:
        scores = [0] * len(mb)
        # pos: position
        # aa: alternate alignment
        for pos, aa in enumerate(mb):
            for variant in aa:
                if variant.is_insertion() or variant.is_deletion():
                    scores[pos] += 2
                elif variant.is_transition():
                    scores[pos] += 1
        if max(scores) > 0:
            # lsi: indices of lower scoring alternate alignments
            lsi = list(x for x in range(len(mb)) if scores[x] < max(scores))
            # remove low scoring alignments from the mismatch block
            # in reverse order so as to not mess up the preceding indices
            lsi.sort(reverse=True)
            for i in lsi:
                mb.pop(i)

    return mb


def prefer_end_of_sorted_list(mb):
    """When all eles fails, prefer the end of the sorted list."""
    for aa in mb:
        aa.sort()
    mb.sort()
    return [mb[-1]]
