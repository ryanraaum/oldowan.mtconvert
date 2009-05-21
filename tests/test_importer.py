from oldowan.mtconvert.importer import load_csv

import os

TEST_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'importer_test_files')
TESTFILE1 = os.path.join(TEST_DIR, 'AbuAmero_et_al_2007.csv')
TESTFILE2 = os.path.join(TEST_DIR, 'Castri_et_al_2009_a.csv')
TESTFILE3 = os.path.join(TEST_DIR, 'Cerny_et_al_2008.csv')
TESTFILE4 = os.path.join(TEST_DIR, 'Non_in_review.csv')
TESTFILE5 = os.path.join(TEST_DIR, 'Richards_et_al_2000_Yemen.csv')

def test_file1():
    popset = load_csv(TESTFILE1,
                      header         = 2,     # number of rows to skip for header info
                      hvr1           = 4,     # column number if present
                      hvr1_covers    = [16000,16569],
                      add16k         = True,  # add 16000 to every hvr1 site?
                      hvr2           = 5,     # column number if present
                      hvr2_covers    = [1,400],
                      sites          = range(6,33), # column number(s) if present
                      sites_on_rCRS  = [3834,6386,6962,7028,8618,8860,8701,10398,10400,10873,11914,11929,12308,12372,12705,14766,15849,15850,15884,15896,15907,15924,15928,15940,15954,15968,15992], # matched entry or list for sites columns
                      haplogroup     = 1,     # column number if present
                      sample_id      = 2,     # column number if present
                      n              = 33,    # column number if present
                      )
                      

def test_file2():
    popset = load_csv(TESTFILE2,
                      header         = 1,     # number of rows to skip for header info
                      hvr1           = 2,     # column number if present
                      hvr1_covers    = [16000,16400],
                      add16k         = True,  # add 16000 to every hvr1 site?
                      rflps          = 3,     # column number if present
                      rflp_format    = 1,     # 1. +/- POSITION ENZYME (i.e. "+16389 HinfI")
                      haplogroup     = 4,     # column number if present
                      sample_id      = 1,     # column number if present
                      haplotype_id   = False, # column number if present
                      n              = 5,     # column number if present
                      population     = 6,     # column number(s) if present
                      )
                      

def test_file3():
    popset = load_csv(TESTFILE3,
                      header         = 1,     # number of rows to skip for header info
                      hvr1           = 2,     # column number if present
                      hvr1_covers    = [16030,16370],
                      add16k         = True,  # add 16000 to every hvr1 site?
                      rflps          = [3,4,5,6,7,8,9], # column number(s) if present
                      rflp_format    = 2,     # 2. '+' present; '-' absent; '#' or '' not tested
                      sites          = False, # column number(s) if present
                      sites_on_rCRS  = False, # matched entry or list for sites columns
                      haplogroup     = 15,    # column number if present
                      sample_id      = False, # column number if present
                      haplotype_id   = 1,     # column number if present
                      pop_with_n     = True,  # are N's arranged by population?
                      n              = [10,11,12,13], # column number(s) if present
                      population     = ["YTA","YTI","YHG","YHA"], # column number or name or if pop_with_n is True, 
                                              # names to go with N columns 
                      )

def test_file4():
    popset = load_csv(TESTFILE4,
                      header         = 1,     # number of rows to skip for header info
                      hvr1           = 2,     # column number if present
                      hvr1_covers    = [16024,16383],
                      add16k         = True,  # add 16000 to every hvr1 site?
                      rflps          = [6,7,8,9,10,11,12,13], # column number if present
                      rflp_format    = 2,     # 2. '+' present; '-' absent; '#' or '' not tested
                      sites          = [14,15], # column number(s) if present
                      sites_on_rCRS  = [769,1018], # matched entry or list for sites columns
                      haplogroup     = 5,     # column number if present
                      sample_id      = 1,     # column number if present
                      sample_id_sep  = ',',   # what separates multiple ids?
                      haplotype_id   = False, # column number if present
                      pop_with_n     = False, # are N's arranged by population?
                      n              = 3,     # column number(s) if present
                      population     = 14,    # column number or name or if pop_with_n is True, 
                                              # names to go with N columns 
                      )

def test_file5():
    popset = load_csv(TESTFILE5,
                      header         = 1,     # number of rows to skip for header info
                      hvr1           = 2,     # column number if present
                      hvr1_covers    = [16090,16365],
                      add16k         = True,  # add 16000 to every hvr1 site?
                      haplogroup     = 1,     # column number if present
                      n              = 3,     # column number(s) if present
                      population     = "Yemenite Jews", # column number or name or if pop_with_n is True, 
                                              # names to go with N columns 
                      )

