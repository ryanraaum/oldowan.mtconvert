from oldowan.mitomotifs import seq2sites
from oldowan.mitomotifs import sites2seq
from oldowan.mitomotifs import str2sites
from oldowan.mitomotifs.polymorphism import Polymorphism

def test_haplotype_1810():
    sites = str2sites('16093C 16183d 16193.1C 16193.2C 16249C')
    seq   = sites2seq(sites, region=range(16000,16570))
    rts   = seq2sites(seq) # rts: round trip sites
    print 'EXP: %s' % sites
    print 'OBS: %s' % rts
    assert sites == rts

def test_haplotype_2236():
    sites = str2sites('16126C 16163G 16185.1T 16185.2T 16189d 16294T 16519C')
    seq   = sites2seq(sites, region=range(16000,16570))
    rts   = seq2sites(seq) # rts: round trip sites
    print 'EXP: %s' % sites
    print 'OBS: %s' % rts
    assert sites == rts

def test_haplotype_2911():
    sites = str2sites('16051G 16129C 16182d 16183d 16193.1C 16193.2C 16362C 16519C')
    seq   = sites2seq(sites, region=range(16000,16570))
    rts   = seq2sites(seq) # rts: round trip sites
    print 'EXP: %s' % sites
    print 'OBS: %s' % rts
    assert sites == rts

def test_haplotype_3070():
    sites = str2sites('16093C 16183d 16184d 16191.1T 16191.2T 16270T')
    seq   = sites2seq(sites, region=range(16000,16570))
    rts   = seq2sites(seq) # rts: round trip sites
    print 'EXP: %s' % sites
    print 'OBS: %s' % rts
    assert sites == rts

def test_haplotype_3805():
    sites = str2sites('16183d 16193.1C 16193.2C 16218T 16519C')
    seq   = sites2seq(sites, region=range(16000,16570))
    rts   = seq2sites(seq) # rts: round trip sites
    print 'EXP: %s' % sites
    print 'OBS: %s' % rts
    assert sites == rts

def test_haplotype_4826():
    sites = str2sites('16172C 16183d 16193.1C 16193.2C 16223T 16320T 16519C')
    seq   = sites2seq(sites, region=range(16000,16570))
    rts   = seq2sites(seq) # rts: round trip sites
    print 'EXP: %s' % sites
    print 'OBS: %s' % rts
    assert sites == rts

def test_haplotype_4827():
    sites = str2sites('16172C 16183d 16193.1C 16193.2C 16223T 16320T')
    seq   = sites2seq(sites, region=range(16000,16570))
    rts   = seq2sites(seq) # rts: round trip sites
    print 'EXP: %s' % sites
    print 'OBS: %s' % rts
    assert sites == rts

