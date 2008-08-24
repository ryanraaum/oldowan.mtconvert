from oldowan.mitomotifs.seq2sites import seq2sites
from oldowan.mitomotifs.polymorphism import Polymorphism

def test_example_1():
    a = Polymorphism(499,0,'A')
    seq = 'ATACTACTAATCTCATCAATACAACCCCCACCCATCCTACCCAGCACACACACACCGCTG'
    result = seq2sites(seq)
    assert len(result) == 1
    assert a in result

def test_example_2():
    a = Polymorphism(498,0,'-')
    b = Polymorphism(499,0,'A')
    seq = 'ATACTACTAATCTCATCAATACAACCCCACCCATCCTACCCAGCACACACACACCGCTG'
    result = seq2sites(seq)
    assert len(result) == 2
    assert a in result
    assert b in result

def test_example_3():
    a = Polymorphism(249,0,'-')
    seq = 'TGCTTGTAGGACATAATAATAACAATTGATGTCTGCACAGCCACTTTCC'
    result = seq2sites(seq)
    assert len(result) == 1
    assert a in result

def test_example_4():
    a = Polymorphism(290,0,'-')
    b = Polymorphism(291,0,'-')
    seq = 'ACACAGACATCATAACAAAATTTCCACCAAACCCCCCC'
    result = seq2sites(seq)
    assert len(result) == 2
    assert a in result
    assert b in result

def test_example_5():
    a = Polymorphism(524,1,'A')
    b = Polymorphism(524,2,'C')
    seq = 'ACAACCCCCGCCCATCCTACCCAGCACACACACACACCGCTGCTAACCCCATACCCC'
    result = seq2sites(seq)
    assert len(result) == 2
    assert a in result
    assert b in result

def test_example_6():
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
    a = Polymorphism(513,0,'-')
    b = Polymorphism(514,0,'-')
    seq = 'ACAACCCCCGCCCATCCTACCCAACACACACACCGCTGCTAACCCCATACCCC'
    result = seq2sites(seq)
    assert len(result) == 2
    assert a in result
    assert b in result

def test_example_8():
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
    a = Polymorphism(16179,0,'T')
    b = Polymorphism(16182,0,'-')
    c = Polymorphism(16183,0,'-')
    d = Polymorphism(16193,1,'C')
    seq = 'CCTGTAGTACATAAAAACCCAATCCACATTAACCCCCTCCCCCATGCTTACAAGCAAGTACAGCAATCAACCCTCAACT'
    result = seq2sites(seq)
    assert len(result) == 4
    assert a in result
    assert b in result
    assert c in result
    assert d in result

def test_example_13():
    a = Polymorphism(16186,0,'T')
    b = Polymorphism(16189,0,'-')
    seq = 'CCTGTAGTACATAAAAACCCAATCCACATCAAAACCTCCCCCCATGCTTACAAGCAAGT'
    result = seq2sites(seq)
    assert len(result) == 2
    assert a in result
    assert b in result

def test_example_14():
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
    a = Polymorphism(16183,0,'-')
    b = Polymorphism(16193,1,'C')
    c = Polymorphism(16193,2,'C')
    seq = 'CCTGTAGTACATAAAAACCCAATCCACATCAAACCCCCTCCCCCCATGCTTACAAGCAAGT'
    result = seq2sites(seq)
    assert len(result) == 3
    assert a in result
    assert b in result
    assert c in result

def test_example_17():
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
    a = Polymorphism(309,0,'-')
    b = Polymorphism(315,1,'C')
    seq = 'CATAACAAAAAATTTCCACCAAACCCCCCTCCCCCCGCTTCTGGCCACAGCACTT'
    result = seq2sites(seq)
    assert len(result) == 2
    assert a in result
    assert b in result

def test_example_19():
    a = Polymorphism(310,1,'T')
    b = Polymorphism(315,1,'C')
    seq = 'AATTTCCACCAAACCCCCCCTTCCCCCCGCTTCTGGCCACAGCACTT'
    result = seq2sites(seq)
    assert len(result) == 2
    assert a in result
    assert b in result

def test_example_20():
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
    assert len(result) == 8
    assert a in result
    assert b in result
    assert c in result
    assert d in result
    assert e in result
    assert f in result
    assert g in result
    assert h in result

