from types import *

from oldowan.mitomotifs.seq2sites import seq2sites
from oldowan.mitomotifs.polymorphism import Polymorphism

from oldowan.mitomotifs.sites2seq import sites2seq

def test_523_524_deletions():
    a = Polymorphism(523,0,'-')
    b = Polymorphism(524,0,'-')
    seq = 'TCTCATCAATACAACCCCCGCCCATCCTACCCAGCACACACACCGCTGCTAACCCCAT'
    result = seq2sites(seq)
    assert a in result
    assert b in result

def test_309_315_insertions():
    a = Polymorphism(309,1,'C')
    b = Polymorphism(315,1,'C')
    seq = 'CATCATAACAAAAAATTTCCACCAAACCCCCCCCTCCCCCCGCTTCTGGCCACAGCACTT'
    result = seq2sites(seq)
    assert a in result
    assert b in result

def test_315_insertion_316_substitution():
    a = Polymorphism(315,1,'C')
    b = Polymorphism(316,0,'A')
    seq = 'CACCAAACCCCCCCTCCCCCCACTTCTGGCCACAGCACTTAAACACATCTCTGCCAAACCCCA'
    result = seq2sites(seq)
    assert a in result
    assert b in result

def test_two_309_insertions():
    a = Polymorphism(309,1,'C')
    b = Polymorphism(309,2,'C')
    seq = 'ATAACAAAAAATTTCCACCAAACCCCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACA' 
    result = seq2sites(seq)
    assert a in result
    assert b in result

def test_two_309_insertions_and_315_insertion():
    a = Polymorphism(309,1,'C')
    b = Polymorphism(309,2,'C')
    c = Polymorphism(315,1,'C')
    seq = 'AATTTCCACCAAACCCCCCCCCTCCCCCCGCTTCTGGCCACAGCACTTAAA'
    result = seq2sites(seq)
    assert len(result) == 3
    assert a in result
    assert b in result
    assert c in result

def test_1378_insertion():
    a = Polymorphism(1378,1,'C')
    seq = 'GGCAAGAAATGGGCTACATTTTCTACCCCCAGAAAACTACGATAGCCCTTATGAAACTTAAGGG'
    result = seq2sites(seq)
    assert a in result

def test_9_insertions():
    seq = 'GCCCGTATTTACCCTATAGCACCCCCTCTAGAGCCCACTGTAAAGCTAACTTAGCATTAACCTT'
    result = seq2sites(seq)
    assert len(result) == 9

def test_498_deletion():
    a = Polymorphism(498,0,'-')
    seq = 'CCCATACTACTAATCTCATCAATACAACCCCGCCCATCCTACCCAG'
    result = seq2sites(seq)
    assert len(result) == 1
    assert a in result

def test_insertion_at_beginning():
    a = Polymorphism(10,0,'G')
    seq = 'GATCACAGGGCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCAT'
    result = seq2sites(seq)
    assert len(result) == 1
    assert a in result

def test_15511_substitution():
    a = Polymorphism(15511,0,'C')
    seq = 'TCACCAGACCTCCTAGGCGACCCAGACAACTATACCCTAGCCAACCCCTTAAA'
    result = seq2sites(seq)
    assert len(result) == 1
    assert a in result

def test_12236_substitution():
    a = Polymorphism(12236,0,'A')
    seq = 'AAAGCTCACAAGAACTGCTAACTCATACCCCCATGTCTAACAACATGGCT'
    result = seq2sites(seq)
    assert len(result) == 1
    assert a in result

def test_3_changes():
    seq = 'GTACATAAAAACCCAATCCACATCAAACCTCCCCCCCATGCTTACAAGCAAGTA'
    result = seq2sites(seq)
    assert len(result) == 3

def test_substitution_at_end():
    a = Polymorphism(16566,0,'A')
    seq = 'CATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACATCACAATG'
    result = seq2sites(seq)
    assert len(result) == 1
    assert a in result

