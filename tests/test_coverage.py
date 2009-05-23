from oldowan.mtconvert.coverage import coverage
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
        seq = sites2seq('', region=item)
        covers = coverage(seq)
        print covers.segments[0].start, item[0]
        print covers.segments[0].stop, item[1]
        assert covers.segments[0].start == item[0]
        assert covers.segments[0].stop  == item[1]

def test_gapped():
    hvr1 = 'TTCTTTCATGGGGAAGCAGATTTGGGTACCACCCAAGTATTGACTCACCCATCAACAACCGCTATGTATTTCGTACATTACTGCCAGCCACCATGAATATTGTACGGTACCATAAATACTTGACCACCTGTAGTACATAAAAACCCAATCCACATCAAAACCCCCTCCCCATGCTTACAAGCAAGTACAGCAATCAACCGTCAACTATCACACATCAACTGCAACTCCAAAGCCACCCCTCACCCACTAGGATACCAACAAACCTACCCACCCTTAACAGTACATAGTACATAAAGCCATTTACCGTACATAGCACATTACAGTCAAATCCCTTCTCGTCCC'
    hvr2 = 'ATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACTTACTAAAGTGTGTTTATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACAC'
    result = coverage(hvr1)
    print result
    result = coverage(hvr2)
    print result
    result = coverage(hvr1+hvr2)
    print result
    #assert False
