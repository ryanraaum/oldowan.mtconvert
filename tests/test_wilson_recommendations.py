from oldowan.mitomotifs.seq2sites import seq2sites
from oldowan.mitomotifs.polymorphism import Polymorphism

def test_example_1():
    """Wilson et al 2002 Example 1

    Seq:  ATACAACCCCCACCCAT

    Seq:  ATACAACCCCCACCCAT
    rCRS: ATACAACCCCCGCCCAT

    Sites: 499A
    """
    a = Polymorphism(499,0,'A')
    seq = 'ATACTACTAATCTCATCAATACAACCCCCACCCATCCTACCCAGCACACACACACCGCTG'
    result = seq2sites(seq)
    assert len(result) == 1
    assert a in result

def test_example_2():
    """Wilson et al 2002 Example 2

    Seq:  ATACAACCCCACCCAT

    Seq:  ATACAACCCC-ACCCAT
    rCRS: ATACAACCCCCGCCCAT

    Sites: 489d 499A
    """
    a = Polymorphism(498,0,'-')
    b = Polymorphism(499,0,'A')
    seq = 'ATACTACTAATCTCATCAATACAACCCCACCCATCCTACCCAGCACACACACACCGCTG'
    result = seq2sites(seq)
    assert len(result) == 2
    assert a in result
    assert b in result

def test_example_3():
    """Wilson et al 2002 Example 3

    Seq:  ATTGATGTC

    Seq:  ATTGA-TGTC
    rCRS: ATTGAATGTC

    Sites: 249d
    """
    a = Polymorphism(249,0,'-')
    seq = 'TGCTTGTAGGACATAATAATAACAATTGATGTCTGCACAGCCACTTTCC'
    result = seq2sites(seq)
    assert len(result) == 1
    assert a in result

def test_example_4():
    """Wilson et al 2002 Example 4

    Seq:  CATAACAAAATTT

    Seq:  CATAACAAAA--TTT
    rCRS: CATAACAAAAAATTT

    Sites: 290d 291d
    """
    a = Polymorphism(290,0,'-')
    b = Polymorphism(291,0,'-')
    seq = 'ACACAGACATCATAACAAAATTTCCACCAAACCCCCCC'
    result = seq2sites(seq)
    assert len(result) == 2
    assert a in result
    assert b in result

def test_example_5():
    """Wilson et al 2002 Example 5

    Seq:  ACCCAGCACACACACACACCGCTG

    Seq:  ACCCAGCACACACACACACCGCTG
    rCRS: ACCCAGCACACACACAC--CGCTG

    Sites: 524.1A 524.2C
    """
    a = Polymorphism(524,1,'A')
    b = Polymorphism(524,2,'C')
    seq = 'ACAACCCCCGCCCATCCTACCCAGCACACACACACACCGCTGCTAACCCCATACCCC'
    result = seq2sites(seq)
    assert len(result) == 2
    assert a in result
    assert b in result

def test_example_6():
    """Wilson et al 2002 Example 6

    Seq:  ACCCAGCACACACCGCTGC

    Seq:  ACCCAGCACACAC----CGCTGC
    rCRS: ACCCAGCACACACACACCGCTGC

    Sites: 521d 522d 523d 524d
    """
    a = Polymorphism(521,0,'-')
    b = Polymorphism(522,0,'-')
    c = Polymorphism(523,0,'-')
    d = Polymorphism(524,0,'-')
    seq = 'ACAACCCCCGCCCATCCTACCCAGCACACACCGCTGCTAACCCCATACCCC'
    result = seq2sites(seq)
    assert len(result) == 4
    assert a in result
    assert b in result
    assert c in result
    assert d in result

def test_example_7():
    """Wilson et al 2002 Example 7

    Seq:  ACCCAACACACACACCGC

    Seq:  ACCCA--ACACACACACCGC
    rCRS: ACCCAGCACACACACACCGC

    Sites: 513d 514d
    """
    a = Polymorphism(513,0,'-')
    b = Polymorphism(514,0,'-')
    seq = 'ACAACCCCCGCCCATCCTACCCAACACACACACCGCTGCTAACCCCATACCCC'
    result = seq2sites(seq)
    assert len(result) == 2
    assert a in result
    assert b in result

def test_example_8():
    """Wilson et al 2002 Example 8

    Seq:  ACCCAGTACACACACCG

    Seq:  ACCCAGTACACACAC--CG
    rCRS: ACCCAGCACACACACACCG

    Sites: 514T 523d 524d
    """
    a = Polymorphism(514,0,'T')
    b = Polymorphism(523,0,'-')
    c = Polymorphism(524,0,'-')
    seq = 'ACAACCCCCGCCCATCCTACCCAGTACACACACCGCTGCTAACCCCATACCCC'
    result = seq2sites(seq)
    assert len(result) == 3
    assert a in result
    assert b in result
    assert c in result

def test_example_9():
    """Wilson et al 2002 Example 9

    Seq:  AAACCCCCCCCTCCCATGCT

    Seq:  AAACCCCCCCCTCCCATGCT
    rCRS: AAAACCCCCTC-CCCATGCT

    Sites: 16183C 15189C 16190.1T
    """
    a = Polymorphism(16183,0,'C')
    b = Polymorphism(16189,0,'C')
    c = Polymorphism(16190,1,'T')
    seq = 'CCTGTAGTACATAAAAACCCAATCCACATCAAACCCCCCCCTCCCATGCTTACAAGCAAGT'
    result = seq2sites(seq)
    assert len(result) == 3
    assert a in result
    assert b in result
    assert c in result

def test_example_10():
    """Wilson et al 2002 Example 10

    Seq:  AACCCCCCCCCCCCATGCT

    Seq:  AACCCCCCCCCCCCATGCT
    rCRS: AAAACCCCCTCCCCATGCT

    Sites: 16182C 16183C 16189C
    """
    a = Polymorphism(16182,0,'C')
    b = Polymorphism(16183,0,'C')
    c = Polymorphism(16189,0,'C')
    seq = 'CCTGTAGTACATAAAAACCCAATCCACATCAACCCCCCCCCCCCATGCTTACAAGCAAGT'
    result = seq2sites(seq)
    assert len(result) == 3
    assert a in result
    assert b in result
    assert c in result

def test_example_11():
    """Wilson et al 2002 Example 11

    Seq:  AAACCCCCCCCCCCCATGCT

    Seq:  AAACCCCCCCCCCCCATGCT
    rCRS: AAAACCCCCTCCCC-ATGCT

    Sites: 16183C 16189C 16193.1C
    """
    a = Polymorphism(16183,0,'C')
    b = Polymorphism(16189,0,'C')
    c = Polymorphism(16193,1,'C')
    seq = 'CCTGTAGTACATAAAAACCCAATCCACATCAAACCCCCCCCCCCCATGCTTACAAGCAAGT'
    result = seq2sites(seq)
    assert len(result) == 3
    assert a in result
    assert b in result
    assert c in result

def test_example_12():
    """Wilson et al 2002 Example 12

    Seq:  TTAACCCCCTCCCCCATGCT

    Seq:  TTAA--CCCCCTCCCCCATGCT
    rCRS: TCAAAACCCCCTCCCC-ATGCT

    Sites: 16179T 16182d 16183d 16193.1C
    """
    a = Polymorphism(16179,0,'T')
    b = Polymorphism(16182,0,'-')
    c = Polymorphism(16183,0,'-')
    d = Polymorphism(16193,1,'C')
    seq = 'CCTGTAGTACATAAAAACCCAATCCACATTAACCCCCTCCCCCATGCTTACAAGCAAGTACAGCAATCAACCCTCAACT'
    result = seq2sites(seq)
    print 'expected: %s' % [a,b,c,d]
    print 'actual:   %s' % result
    assert len(result) == 4
    assert a in result
    assert b in result
    assert c in result
    assert d in result

def test_example_13():
    """Wilson et al 2002 Example 13

    Seq:  AAAACCTCCCCCCATGCT

    Seq:  AAAACCTCC-CCCCATGCT
    rCRS: AAAACCCCCTCCCCATGCT

    Sites: 16186T 16189d
    """
    a = Polymorphism(16186,0,'T')
    b = Polymorphism(16189,0,'-')
    seq = 'CCTGTAGTACATAAAAACCCAATCCACATCAAAACCTCCCCCCATGCTTACAAGCAAGT'
    result = seq2sites(seq)
    assert len(result) == 2
    assert a in result
    assert b in result

def test_example_14():
    """Wilson et al 2002 Example 14

    Seq:  AAACCCCCCCTCCCCCATGCT

    Seq:  AAACCCCCCCTCCCCCATGCT
    rCRS: AAAACCCCC-TCCCC-ATGCT

    Sites: 16183C 16188.1C 16193.1C
    """
    a = Polymorphism(16183,0,'C')
    b = Polymorphism(16188,1,'C')
    c = Polymorphism(16193,1,'C')
    seq = 'CCTGTAGTACATAAAAACCCAATCCACATCAAACCCCCCCTCCCCCATGCTTACAAGCAAGT'
    result = seq2sites(seq)
    assert len(result) == 3
    assert a in result
    assert b in result
    assert c in result

def test_example_15():
    """Wilson et al 2002 Example 15

    Seq:  TTAAACCCCCCCCTCCCATGCT

    Seq:  TTAAACCCCCCCCTCCCATGCT
    rCRS: TCAAAACCCCCTC-CCCATGCT

    Sites: 16179T 16183C 16189C 16190.1T
    """
    a = Polymorphism(16179,0,'T')
    b = Polymorphism(16183,0,'C')
    c = Polymorphism(16189,0,'C')
    d = Polymorphism(16190,1,'T')
    seq = 'CCTGTAGTACATAAAAACCCAATCCACATTAAACCCCCCCCTCCCATGCTTACAAGCAAGTACAGCAATCAACCCTCAACT'
    result = seq2sites(seq)
    assert len(result) == 4
    assert a in result
    assert b in result
    assert c in result
    assert d in result

def test_example_16():
    """Wilson et al 2002 Example 16

    Seq:  AAACCCCCTCCCCCCATGCT

    Seq:  AAA-CCCCCTCCCCCCATGCT
    rCRS: AAAACCCCCTCCCC--ATGCT

    Sites: 16183d 16193.1C 16193.2C
    """
    a = Polymorphism(16183,0,'-')
    b = Polymorphism(16193,1,'C')
    c = Polymorphism(16193,2,'C')
    seq = 'CCTGTAGTACATAAAAACCCAATCCACATCAAACCCCCTCCCCCCATGCTTACAAGCAAGT'
    result = seq2sites(seq)
    print 'expected: %s' % [a,b,c]
    print 'actual:   %s' % result
    assert len(result) == 3
    assert a in result
    assert b in result
    assert c in result
 
def test_example_17():
    """Wilson et al 2002 Example 17

    Seq:  AAACCCCCCCCCGC

    Seq:  AAACCCCCCC----CCGC
    rCRS: AAACCCCCCCTCCCCCGC

    Sites: 310d 311d 312d 313d
    """
    a = Polymorphism(310,0,'-')
    b = Polymorphism(311,0,'-')
    c = Polymorphism(312,0,'-')
    d = Polymorphism(313,0,'-')
    seq = 'AATTTCCACCAAACCCCCCCCCGCTTCTGGCCACAGCACTT'
    result = seq2sites(seq)
    assert len(result) == 4
    assert a in result
    assert b in result
    assert c in result
    assert d in result

def test_example_18():
    """Wilson et al 2002 Example 18

    Seq:  AAACCCCCCTCCCCCCGC

    Seq:  AAACCCCCC-TCCCCCCGC
    rCRS: AAACCCCCCCTCCCCC-GC

    Sites: 309d 315.1C
    """
    a = Polymorphism(309,0,'-')
    b = Polymorphism(315,1,'C')
    seq = 'CATAACAAAAAATTTCCACCAAACCCCCCTCCCCCCGCTTCTGGCCACAGCACTT'
    result = seq2sites(seq)
    print 'expected: %s' % [a,b]
    print 'actual:   %s' % result
    assert len(result) == 2
    assert a in result
    assert b in result

def test_example_19():
    """Wilson et al 2002 Example 19

    Seq:  AAACCCCCCCTTCCCCCCGC

    Seq:  AAACCCCCCCTTCCCCCCGC
    rCRS: AAACCCCCCCT-CCCCC-GC

    Sites: 310.1T 315.1C
    """
    a = Polymorphism(310,1,'T')
    b = Polymorphism(315,1,'C')
    seq = 'AATTTCCACCAAACCCCCCCTTCCCCCCGCTTCTGGCCACAGCACTT'
    result = seq2sites(seq)
    assert len(result) == 2
    assert a in result
    assert b in result

def test_example_20():
    """Wilson et al 2002 Example 20

    Seq:  AAAGACACCCCCCCCCCCCACA

    Seq:  AAAGACACCCCCCCCCCCCACA
    rCRS: AAAGACACCCCCC------ACA

    Sites: 573.1C 573.2C 573.3C 573.4C 573.5C 573.6C
    """
    a = Polymorphism(573,1,'C')
    b = Polymorphism(573,2,'C')
    c = Polymorphism(573,3,'C')
    d = Polymorphism(573,4,'C')
    e = Polymorphism(573,5,'C')
    f = Polymorphism(573,6,'C')
    seq = 'AACCAAACCCCAAAGACACCCCCCCCCCCCACAGTTTATGTAGCTT'
    result = seq2sites(seq)
    assert len(result) == 6
    assert a in result
    assert b in result
    assert c in result
    assert d in result
    assert e in result
    assert f in result

def test_example_21():
    """Wilson et al 2002 Example 21

    Seq:  CTGGAGCACCC

    Seq:  CTGGAGC------ACCC
    rCRS: CTGGAGCCGGAGCACCC

    Sites: 105d 106d 107d 108d 109d 110d
    """
    a = Polymorphism(105,0,'-')
    b = Polymorphism(106,0,'-')
    c = Polymorphism(107,0,'-')
    d = Polymorphism(108,0,'-')
    e = Polymorphism(109,0,'-')
    f = Polymorphism(110,0,'-')
    seq = 'GCATTGCGAGACGCTGGAGCACCCTATGTCGCAGTATCT'
    result = seq2sites(seq)
    assert len(result) == 6
    assert a in result
    assert b in result
    assert c in result
    assert d in result
    assert e in result
    assert f in result

def test_example_22():
    """Wilson et al 2002 Example 22

    Seq:  AGATCCTGGAGCCCCC

    Seq:  AGATC-CTGGAGCC------CCC
    rCRS: AGA-CGCTGGAGCCGGAGCACCC

    Sites: 95.1T 97d 106d 107d 108d 109d 110d 111d
    """
    a = Polymorphism(95,1,'T')
    b = Polymorphism(97,0,'-')
    c = Polymorphism(106,0,'-')
    d = Polymorphism(107,0,'-')
    e = Polymorphism(108,0,'-')
    f = Polymorphism(109,0,'-')
    g = Polymorphism(110,0,'-')
    h = Polymorphism(111,0,'-')
    seq = 'TCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGATCCTGGAGCCCCCTATGTCGCAGTATCT'
    result = seq2sites(seq)
    print 'expected: %s' % [a,b,c,d,e,f,g,h]
    print 'actual:   %s' % result
    assert len(result) == 8
    assert a in result
    assert b in result
    assert c in result
    assert d in result
    assert e in result
    assert f in result
    assert g in result
    assert h in result

