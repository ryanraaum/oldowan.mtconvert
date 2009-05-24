from oldowan.mtconvert.coverage import calc_coverage
from oldowan.mtconvert.coverage import Coverage
from oldowan.mtconvert.sites2seq import sites2seq

def test_single_rcrs_segment():
    ranges = [(16024,16365),
              (73,340),
              (16000,340),
              (16000,577),
              (577,15992),
              (16000,16569),
             ]
    for item in ranges:
        should_cover = Coverage(item)
        seq = sites2seq('', region=should_cover)
        covers = calc_coverage(seq)
        print "should cover: '%s' - does cover: '%s'" % (should_cover, covers)
        assert should_cover == covers

def test_gapped_rcrs_segments():
    ranges = [ [(16024,16365), (73,340) ],
               [(16024,16569), (73,340) ],
               [(16024,16569), (73,577) ],
               [(16024,16400), (1,577) ],
             ]
    for pair in ranges:
        should_cover = Coverage(*pair)
        seq = sites2seq('', region=should_cover)
        covers = calc_coverage(seq)
        print "should cover: '%s' - does cover: '%s'" % (should_cover, covers)
        assert should_cover == covers

def test_single_variant_segments():
    should_cover = Coverage((16024,16365))
    variants = ['16311',
                '16223 16265T',
                '16223 16260 16270 16287 16311',
                '16148 16172 16187 16188G 16189 16209 16223 16230 16311 16320',
                '16024', # tricky, first site is a variant
                '16365', # tricky, last site is a variant
                '16024 16365', # tricky, both first and last sites are variants
                '16363 16364 16365', # tricky, last few sites are variants
               ]
    for v in variants:
        seq = sites2seq(v, region=should_cover)
        covers = calc_coverage(seq)
        print "should cover: '%s' - does cover: '%s'" % (should_cover, covers)
        assert should_cover == covers

def test_gapped_variant_segments():
    should_cover = Coverage((16024,16365),(73,340))
    variants = ['16311 263 311',
                '16024 73', # little bit tricky - first site from second chunk is variant 
                # '16365', # impossible, middle site is a variant
               ]
    for v in variants:
        seq = sites2seq(v, region=should_cover)
        covers = calc_coverage(seq)
        print "should cover: '%s' - does cover: '%s'" % (should_cover, covers)
        assert should_cover == covers

def test_more_than_two_gapped_segments():
    sets = [ [(16024,16365),(73,340),(577,1000)],
             [(16024,16365),(16300,340),(577,1000)], # overlapping segments
           ]
    for s in sets:
        should_cover = Coverage(*s)
        seq = sites2seq('', region=should_cover)
        covers = calc_coverage(seq)
        print "should cover: '%s' - does cover: '%s'" % (should_cover, covers)
        assert should_cover == covers


