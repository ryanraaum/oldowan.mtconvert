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
        print should_cover, covers
        assert should_cover == covers

