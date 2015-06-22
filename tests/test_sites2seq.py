from oldowan.mtconvert.sites2seq import sites2seq
from oldowan.mtdna import rCRS
from oldowan.polymorphism import Polymorphism

def test_single_site_as_string():
    sites = '16129-'
    seq   = sites2seq(sites)
    assert '-' in seq

def test_single_site_as_string_in_list():
    sites = ['16129-']
    seq   = sites2seq(sites)
    assert '-' in seq

def test_single_site_as_Polymorphism_in_list():
    sites = [Polymorphism(16129,0,'-')]
    seq   = sites2seq(sites)
    assert '-' in seq

def test_two_sites_in_string():
    sites = '316C 316.1A'
    seq = sites2seq(sites, region=range(314,319))
    assert seq == 'CCCACT'

def test_two_sites_as_string_in_list():
    sites = ['316C', '316.1A']
    seq = sites2seq(sites, region=(314,318))
    assert seq == 'CCCACT'

def test_two_sites_as_Polymorphisms_in_list():
    sites = [Polymorphism(316,0,'C'), Polymorphism(316,1,'A')]
    seq = sites2seq(sites, region=(314,318))
    assert seq == 'CCCACT'

def test_mix_of_string_and_Polymorphisms_in_list():
    sites = [Polymorphism(316,0,'C'), '316.1A']
    seq = sites2seq(sites, region=(314,318))
    assert seq == 'CCCACT'
    # reverse the order
    sites = ['316.1A', Polymorphism(316,0,'C')]
    seq = sites2seq(sites, region=(314,318))
    assert seq == 'CCCACT'

def test_default_sequence_range_is_HVR1():
    assert sites2seq('') == rCRS[16023:16365]
    assert sites2seq('', region=(16024,16365)) == rCRS[16023:16365]

def test_rCRS_argument_returns_rCRS():
    assert sites2seq('rCRS') == rCRS[16023:16365]

def test_HVR2_range_argument():
    assert sites2seq('', region='HVR2') == rCRS[72:340]

def test_HVR1andHVR2_range_argument():
    assert sites2seq('', region='HVR1and2') == rCRS[16023:16365]+rCRS[72:340]

def test_HVR1toHVR2_range_argument():
    assert sites2seq('', region='HVR1to2') == rCRS[16023:]+rCRS[:340]

def test_coding_range_argument():
    assert sites2seq('', region='coding') == rCRS[576:15992]

def test_all_range_argument():
    assert sites2seq('', region='all') == rCRS

def test_single_substition_in_correct_place():
    # rCRS string is 0-based count
    # rCRS as DNA sequence is 1-based
    # the the numbering discrepancy
    assert 'A' == sites2seq('16129A', region='all')[16128]

def test_custom_range():
    assert 'GAG' == sites2seq('', region=[1,5,9])

def test_insertion_in_custom_range():
    assert 'GCAT' == sites2seq('1.1C', region=[1,2,3])

def test_two_insertions_in_custom_range():
    assert 'GCAGT' == sites2seq('1.1C 2.1G', region=[1,2,3])

def test_deletion_at_start_of_custom_range():
    assert '-AT' == sites2seq('1-', region=[1,2,3])

def test_deletion_in_middle_of_custom_range():
    assert 'G-T' == sites2seq('2-', region=[1,2,3])

def test_sites_outside_region_snp():
    assert sites2seq('73G') == rCRS[16023:16365]

def test_sites_outside_region_del():
    assert sites2seq('489d') == rCRS[16023:16365]

def test_sites_outside_region_ins():
    assert sites2seq('315.1C') == rCRS[16023:16365]
