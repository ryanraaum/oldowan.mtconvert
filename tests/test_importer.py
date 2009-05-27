from oldowan.mtconvert.importer import load_csv
from oldowan.mtconvert.coverage import Coverage

import os

TEST_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'importer_test_files')
TESTFILE1 = os.path.join(TEST_DIR, 'AbuAmero_et_al_2007.csv')
TESTFILE2 = os.path.join(TEST_DIR, 'Castri_et_al_2009.csv')
TESTFILE3 = os.path.join(TEST_DIR, 'Cerny_et_al_2008.csv')
TESTFILE4 = os.path.join(TEST_DIR, 'Non_in_review.csv')
TESTFILE5 = os.path.join(TEST_DIR, 'Richards_et_al_2000_Yemen.csv')

def test_file1():
    popset = load_csv(TESTFILE1,
                      header         = 2,     # number of rows to skip for header info
                      hvr1           = 3,     # column number if present
                      hvr1_covers    = [16000,16569],
                      add16k         = True,  # add 16000 to every hvr1 site?
                      hvr2           = 4,     # column number if present
                      hvr2_covers    = [1,400],
                      sites          = range(5,32), # column number(s) if present
                      sites_on_rCRS  = [3834,6386,6962,7028,8618,8860,8701,10398,10400,10873,11914,11929,12308,12372,12705,14766,15849,15850,15884,15896,15907,15924,15928,15940,15954,15968,15992], # matched entry or list for sites columns
                      haplogroup     = 0,     # column number if present
                      sample_id      = 1,     # column number if present
                      n              = 32,    # column number if present
                      population     = 'Saudi Arabia',
                      )

    for line, error in popset.errors:
        print "%d: %s: %s" % (line, error.message, error.expression)
    assert not popset.errors

    # the coverage for the population is the intersection of all the coverages of all the
    # samples in that population. In this case, none of the additional sites were typed
    # in ALL of the samples
    should_cover=Coverage((16000,16569),(1,400))
    print "should cover: '%s' - does cover: '%s'" % (should_cover, popset.coverage)
    assert popset.coverage                  == should_cover
    assert popset.populations[0].coverage   == should_cover

    # make sure that a selection of the samples have their appropriate coverage
    # sample '27' was typed for 7028, 12308, 12372 as well
    should_cover=Coverage((16000,16569),(1,400),7028,12308,12372)
    sample_coverage = popset.populations[0].sample_by_id('27').coverage
    print "should cover: '%s' - does cover: '%s'" % (should_cover, sample_coverage)
    assert sample_coverage == should_cover

    # sample '223' was typed for only 7028 as well
    should_cover=Coverage((16000,16569),(1,400),7028)
    sample_coverage = popset.populations[0].sample_by_id('223').coverage
    print "should cover: '%s' - does cover: '%s'" % (should_cover, sample_coverage)
    assert sample_coverage == should_cover

    # sample '50' was typed for 11929, 12308, 12372 as well
    should_cover=Coverage((16000,16569),(1,400),11929,12308,12372)
    sample_coverage = popset.populations[0].sample_by_id('50').coverage
    print "should cover: '%s' - does cover: '%s'" % (should_cover, sample_coverage)
    assert sample_coverage == should_cover

    # sample '201' was typed for 6962, 7028, 10398, 10400, 10873, 12705, 15884, 15896 as well
    should_cover=Coverage((16000,16569),(1,400),6962,7028,10398,10400,10873,12705,15884,15896)
    sample_coverage = popset.populations[0].sample_by_id('201').coverage
    print "should cover: '%s' - does cover: '%s'" % (should_cover, sample_coverage)
    assert sample_coverage == should_cover

    assert popset.num_populations == 1
    assert popset.num_samples     == 120
                      
    pop = popset.populations[0]
    haps = ['H'] * 15 + ['I'] + ['J'] * 7 + ['J1'] * 4 + ['J1b'] * 14 + ['J1d'] * 2 + \
    ['K'] * 7 + ['L2a1a','L2a2','L2c2','L3d1','L3f','L3h1','L3h1','L3i','M'] + \
    ['M1a'] * 4 + ['M1b1','M25','M3'] + ['N1a'] * 5 + ['N1b'] * 3 + ['N1c'] * 2 + \
    ['preHV'] * 21 + ['T','T1'] + ['T3'] * 4 + \
    ['T5'] * 2 + ['U1a'] * 2 + ['U1b','U2e'] + ['U3'] * 2 + \
    ['U5a1a','U6a','U8b','U9','U9','U9','W','X','X']
    
    for i in range(popset.num_samples):
        print pop.samples[i].haplogroup, haps[i]
        assert pop.samples[i].haplogroup == haps[i]
        # check that the population name is correctly assigned while we're looping 
        # through the samples anyways
        assert pop.samples[i].population == 'Saudi Arabia'


def test_file2():
    popset = load_csv(TESTFILE2,
                      header         = 1,     # number of rows to skip for header info
                      hvr1           = 1,     # column number if present
                      hvr1_covers    = [16000,16400],
                      add16k         = True,  # add 16000 to every hvr1 site?
                      rflps          = 2,     # column number if present
                      rflp_format    = 1,     # 1. +/- POSITION ENZYME (i.e. "+16389 HinfI")
                      haplogroup     = 3,     # column number if present
                      sample_id      = 0,     # column number if present
                      haplotype_id   = False, # column number if present
                      n              = 4,     # column number if present
                      population     = 5,     # column number(s) if present
                      )
                      

def test_file3():
    popset = load_csv(TESTFILE3,
                      header         = 1,     # number of rows to skip for header info
                      hvr1           = 1,     # column number if present
                      hvr1_covers    = [16030,16370],
                      add16k         = True,  # add 16000 to every hvr1 site?
                      rflps          = [2,3,4,5,6,7,8], # column number(s) if present
                      rflp_format    = 2,     # 2. '+' present; '-' absent; '#' or '' not tested
                      sites          = False, # column number(s) if present
                      sites_on_rCRS  = False, # matched entry or list for sites columns
                      haplogroup     = 14,    # column number if present
                      sample_id      = False, # column number if present
                      haplotype_id   = 0,     # column number if present
                      pop_with_n     = True,  # are N's arranged by population?
                      n              = [9,10,11,12], # column number(s) if present
                      population     = ["YTA","YTI","YHG","YHA"], # column number or name or if pop_with_n is True, 
                                              # names to go with N columns 
                      )

def test_file4():
    popset = load_csv(TESTFILE4,
                      header         = 1,     # number of rows to skip for header info
                      hvr1           = 1,     # column number if present
                      hvr1_covers    = [16024,16383],
                      add16k         = True,  # add 16000 to every hvr1 site?
                      rflps          = range(5,13), # column number if present
                      rflp_format    = 2,     # 2. '+' present; '-' absent; '#' or '' not tested
                      sites          = [13,14], # column number(s) if present
                      sites_on_rCRS  = [769,1018], # matched entry or list for sites columns
                      haplogroup     = 4,     # column number if present
                      sample_id      = 0,     # column number if present
                      sample_id_sep  = ',',   # what separates multiple ids?
                      haplotype_id   = False, # column number if present
                      pop_with_n     = False, # are N's arranged by population?
                      n              = 2,     # column number(s) if present
                      population     = 15,    # column number or name or if pop_with_n is True, 
                                              # names to go with N columns 
                      )

def test_file5():
    popset = load_csv(TESTFILE5,
                      header         = 1,     # number of rows to skip for header info
                      hvr1           = 1,     # column number if present
                      hvr1_covers    = [16090,16365],
                      add16k         = True,  # add 16000 to every hvr1 site?
                      haplogroup     = 0,     # column number if present
                      n              = 2,     # column number(s) if present
                      population     = "Yemenite Jews", # column number or name or if pop_with_n is True, 
                                              # names to go with N columns 
                      )

