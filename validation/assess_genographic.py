import os
here = os.path.abspath(os.path.dirname(__file__))
module_dir = os.path.join(here, '..')
genographic_sites_filepath = os.path.join(here, 'genographic_sites.txt')

import sys
sys.path = [module_dir] + sys.path

outfn = os.path.join(here, 'fail_genographic.txt')
outf = open(outfn, 'w')

from oldowan.mitomotifs import sites2seq
from oldowan.mitomotifs import seq2sites
from oldowan.mitomotifs import str2sites

count = 0
for line in open(genographic_sites_filepath, 'rU'):
    count += 1
    line = line.upper()
    sites = str2sites(line)
    sites.sort()
    seq = sites2seq(sites, region=range(16000,16570))
    roundtrip_sites = seq2sites(seq)
    roundtrip_sites.sort()
    if sites != roundtrip_sites:
        outf.write("%d,%s,%s\n" % (count, sites, roundtrip_sites))
    else:
        print "%d,good" % count
    outf.flush()

outf.close()
