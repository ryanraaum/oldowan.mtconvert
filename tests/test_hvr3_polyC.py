from oldowan.mtconvert.seq2sites import seq2sites
from oldowan.polymorphism import Polymorphism

def test_normal_polyC():
    """Normal Poly C stretch at end of HVR3

    Seq:  CAAAGACACCCCCCACA

    Seq:  CAAAGACACCCCCCACA
    rCRS: CAAAGACACCCCCCACA

    Sites: <None>
    """
    seq = 'CAAAGACACCCCCCACA'
    result = seq2sites(seq)
    assert len(result) == 0

def test_expanded_polyC():
    """Expanded Poly C stretch at end of HVR3

    Seq:  CAAAGACACCCCCCCCCACA

    Seq:  CAAAGACACCCCCCCCCACA
    rCRS: CAAAGACACCCCCC---ACA

    Sites: 573.1C 573.2C 573.3C
    """
    a = Polymorphism(573,1,'C')
    b = Polymorphism(573,2,'C')
    c = Polymorphism(573,3,'C')
    d = Polymorphism(573,4,'C')
    e = Polymorphism(573,5,'C')
    f = Polymorphism(573,6,'C')
    seq = 'ACCCCATACCCCGAACCAACCAAACCCCAAAGACACCCCCCCCCCCCACA'
    result = seq2sites(seq)
    assert len(result) == 6
    assert a in result
    assert b in result
    assert c in result
    assert d in result
    assert e in result
    assert f in result
