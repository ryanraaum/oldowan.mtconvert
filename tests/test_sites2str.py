from oldowan.mitomotifs.sites2str import sites2str
from oldowan.mitomotifs.polymorphism import Polymorphism

def test_no_sites():
    assert '' == sites2str([])

def test_unnested():
    a = Polymorphism(16129,0,'A')
    b = Polymorphism(16456,1,'G')
    assert '16129A 16456.1G' == sites2str([a, b])
    
def test_one_unresolved_ambiguity():
    a = Polymorphism(16129,0,'A')
    b = Polymorphism(16455,1,'G')
    c = Polymorphism(16456,1,'G')
    assert '16129A (16455.1G or 16456.1G)' == sites2str([a, [[b], [c]]])

def test_several_unresolved_ambiguities():
    a = Polymorphism(16129,1,'A')
    b = Polymorphism(16130,1,'G')
    c = Polymorphism(16129,1,'A')
    d = Polymorphism(16129,2,'G')
    e = Polymorphism(16128,1,'G')
    f = Polymorphism(16128,2,'A')
    g = Polymorphism(16223,0,'G')
    expected = '((16129.1A 16130.1G) or (16129.1A 16129.2G) or (16128.1G 16128.2A)) 16223G'
    sites = [[[a, b], [c, d], [e, f]], g]
    assert expected == sites2str(sites)

