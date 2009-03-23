from oldowan.polymorphism import Polymorphism

from oldowan.mtconvert.reduction_funcs import prefer_fewer
from oldowan.mtconvert.reduction_funcs import prefer_indels_at_end
from oldowan.mtconvert.reduction_funcs import prefer_multi_inserts

def test_fewer_is_better_one_vs_two():
    input  = [['a'], ['b', 'c']]
    output = [['a']]
    assert output == prefer_fewer(input)

def test_fewer_is_better_zero_vs_one_vs_two():
    input  = [['a'], ['b', 'c'], []]
    output = [[]]
    assert output == prefer_fewer(input)
    
def test_fewer_is_better_two_with_one_and_one_with_two():
    input  = [['a'], ['b', 'c'], ['d']]
    output = [['a'], ['d']]
    assert output == prefer_fewer(input)

def test_two_consecutive_insertions():
    a = Polymorphism(1,1,'A')
    b = Polymorphism(2,1,'A')
    input  = [[a], [b]]
    output = [[b]]
    assert output == prefer_indels_at_end(input)

def test_pick_most_3_prime_insertion():
    # alternates are several contiguous insertions
    a = Polymorphism(1,1,'A')
    b = Polymorphism(2,1,'A')
    c = Polymorphism(3,1,'A')
    d = Polymorphism(4,1,'A')
    e = Polymorphism(5,1,'A')
    input  = [[a], [e], [b], [c], [d]]
    output = [[e]]
    assert output == prefer_indels_at_end(input)

def test_pick_most_3_prime_pair_of_insertions():
    # alternates are two contiguous insertions, with a shared site
    a = Polymorphism(1,1,'A')
    b = Polymorphism(5,0,'A')
    c = Polymorphism(2,1,'A')
    d = Polymorphism(5,0,'A')
    input  = [[a,b], [c,d]]
    output = [[c,d]]
    assert output == prefer_indels_at_end(input)

def test_pick_most_3_prime_pair_of_separated_insertions():
    # alternates are two near, but not contiguous pairs
    a = Polymorphism(1,1,'A')
    b = Polymorphism(2,1,'A')
    c = Polymorphism(2,1,'A')
    d = Polymorphism(3,1,'A')
    input = [[a,b],[c,d]]
    output = [[c,d]]
    result = prefer_indels_at_end(input)
    for x in result:
        x.sort()
    assert output == result

def test_pick_most_3_prime_double_insertion():
    # pairs of multi-inserts
    a = Polymorphism(1,1,'A')
    b = Polymorphism(1,2,'A')
    c = Polymorphism(2,1,'A')
    d = Polymorphism(2,2,'A')
    e = Polymorphism(3,1,'A')
    f = Polymorphism(3,2,'A')
    input = [[a,b],[c,d],[e,f]]
    output = [[e,f]]
    assert output == prefer_indels_at_end(input)

def test_pick_most_3_prime_set_of_insertions():
    # triplets of multi-inserts
    a = Polymorphism(1,1,'A')
    b = Polymorphism(1,2,'A')
    c = Polymorphism(5,1,'A')
    d = Polymorphism(2,1,'A')
    e = Polymorphism(2,2,'A')
    f = Polymorphism(6,1,'A')
    input = [[a,b,c],[d,e,f]]
    output = [[d,e,f]]
    assert output == prefer_indels_at_end(input)

def test_pick_most_3_prime_single_deletion():
    # alternates are two contiguous deletions
    a = Polymorphism(1,0,'-')
    b = Polymorphism(2,0,'-')
    input  = [[a], [b]]
    output = [[b]]
    assert output == prefer_indels_at_end(input)

def test_pick_most_3_prime_of_many_single_deletions():
    # alternates are several contiguous deletions
    a = Polymorphism(1,0,'-')
    b = Polymorphism(2,0,'-')
    c = Polymorphism(3,0,'-')
    d = Polymorphism(4,0,'-')
    e = Polymorphism(5,0,'-')
    input  = [[a], [e], [b], [c], [d]]
    output = [[e]]
    assert output == prefer_indels_at_end(input)

def test_pick_most_3_prime_deletion_carrying_a_substitution():
    # alternates are two contiguous deletions, with a shared site
    a = Polymorphism(1,0,'-')
    b = Polymorphism(2,0,'-')
    c = Polymorphism(5,0,'A')
    d = Polymorphism(5,0,'A')
    input  = [[a,c], [b,d]]
    output = [[b,d]]
    assert output == prefer_indels_at_end(input)

def test_most_3_prime_matched_pair_of_deletions():
    # alternates are 2 pairs of 2 contiguous deletions
    a = Polymorphism(1,0,'-')
    b = Polymorphism(2,0,'-')
    c = Polymorphism(20,0,'-')
    d = Polymorphism(21,0,'-')
    e = Polymorphism(1,0,'-')
    f = Polymorphism(2,0,'-')
    g = Polymorphism(20,0,'-')
    h = Polymorphism(21,0,'-')
    input = [[a,c], [e,d], [b,g], [f,h]]
    output = [[f,h]]
    assert output == prefer_indels_at_end(input)

def test_pick_most_3_prime_pair_of_deletions_from_many():
    # alternates are 2 pairs of contiguous deletions
    # one is a choice among 2, the other is a choice among 3
    a = Polymorphism(1,0,'-')
    b = Polymorphism(20,0,'-')
    c = Polymorphism(1,0,'-')
    d = Polymorphism(21,0,'-')
    e = Polymorphism(2,0,'-')
    f = Polymorphism(20,0,'-')
    g = Polymorphism(2,0,'-')
    h = Polymorphism(21,0,'-')
    i = Polymorphism(3,0,'-')
    j = Polymorphism(20,0,'-')
    k = Polymorphism(3,0,'-')
    l = Polymorphism(21,0,'-')
    input = [[a,b], [c,d], [e,f], [g,h], [i,j], [k,l]]
    output = [[k,l]]
    assert output == prefer_indels_at_end(input)

def test_that_multiple_deletions_in_an_alignment_are_not_lost():
    # two deletions in a row in the same alt_aln should not be reduced
    a = Polymorphism(1,0,'-')
    b = Polymorphism(2,0,'-')
    input = [[a,b]]
    output = [[a,b]]
    assert output == prefer_indels_at_end(input)
    
def test_pick_3_prime_pair_of_deletions():
    # paired deletions should slide to the end
    a = Polymorphism(1,0,'-')
    b = Polymorphism(2,0,'-')
    c = Polymorphism(2,0,'-')
    d = Polymorphism(3,0,'-')
    input = [[a,b],[c,d]]
    output = [[c,d]]
    result = prefer_indels_at_end(input)
    for x in result:
        x.sort()
    assert output == result

def test_pick_3_prime_set_of_deletions():
    # set of deletions should slide to the end
    a = Polymorphism(1,0,'-')
    b = Polymorphism(2,0,'-')
    c = Polymorphism(3,0,'-')
    d = Polymorphism(4,0,'-')
    e = Polymorphism(2,0,'-')
    f = Polymorphism(3,0,'-')
    g = Polymorphism(4,0,'-')
    h = Polymorphism(5,0,'-')
    input = [[a,b,c,d],[e,f,g,h]]
    output = [[e,f,g,h]]
    result = prefer_indels_at_end(input)
    for x in result:
        x.sort()
    assert output == prefer_indels_at_end(input)

def test_pick_3_prime_set_of_deletions_and_substitutions():
    # set of deletions should slide to the end
    a = Polymorphism(27,0,'T')
    b = Polymorphism(29,0,'T')
    c = Polymorphism(33,0,'-')
    d = Polymorphism(30,0,'T')
    e = Polymorphism(33,0,'C')
    f = Polymorphism(27,0,'-')
    input = [[a,b,c],[d,e,f]]
    output = [[a,b,c]]
    result = prefer_indels_at_end(input)
    for x in result:
        x.sort()
    assert output == prefer_indels_at_end(input)

def test_prefer_multi_inserts_does_not_affect_single_insertion_alternatives():
    # make sure non-cases make it through
    a = Polymorphism(1,0,'A')
    b = Polymorphism(1,1,'A')
    c = Polymorphism(1,0,'A')
    d = Polymorphism(2,1,'A')
    input = [[a,b],[c,d]]
    output = [[a,b],[c,d]]
    assert output == prefer_multi_inserts(input)

def test_pick_one_two_base_insertion_over_two_separate_insertions():
    a = Polymorphism(1,1,'A')
    b = Polymorphism(1,2,'A')
    c = Polymorphism(1,1,'A')
    d = Polymorphism(2,1,'A')
    input = [[a,b],[c,d]]
    output = [[a,b]]
    assert output == prefer_multi_inserts(input)
    

