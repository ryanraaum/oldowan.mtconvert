import os
here = os.path.abspath(os.path.dirname(__file__))
module_dir = os.path.join(here, '..')

import sys
sys.path = [module_dir] + sys.path

from oldowan.mitomotifs import seq2sites
from oldowan.mitomotifs import sites2seq

mtgs_fn = os.path.join(here, 'all_mt_genomes.txt')
ofn = os.path.join(here, 'failures.txt')

# not mitochondrial, too many N's, not complete, etc.
expect_to_fail = [393, 414, 1327, 1328, 1329, 1330, 1331, 1332, 1333, 
                  1346, 1347, 2195, 2196, 2199, 2206, 2207, 2209, 2210, 2211, 
                  2212, 2216, 2217, 2218, 2219, 2220, 2224, 2227, 2228, 2229, 
                  2232, 2234, 2236, 2238, 2244, 2278, 2279, 2281, 2282, 2288, 
                  2667, 2668, 3186, 3189, 3916, 3947, 4089, 4092 ]

count = 0
of = open(ofn, 'w')
for entry in open(mtgs_fn, 'U'):
    count += 1
    if count not in expect_to_fail:
        entry = entry.strip()
        entry = entry.replace('-', '')
        try:
            sites = seq2sites(entry)
            seq = sites2seq(sites, what='all')
            seq = seq.replace('-', '')
            if entry != seq:
                print count, 'FAILED'
                of.write('%d\n' % count)
            else:
                print count, 'PASSED'
        except:
            print count, 'FAILED'
            of.write('%d\n' % count)
        of.flush()
of.close()
