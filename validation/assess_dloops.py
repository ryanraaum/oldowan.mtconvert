import os
here = os.path.abspath(os.path.dirname(__file__))
module_dir = os.path.join(here, '..')

import sys
sys.path = [module_dir] + sys.path

from oldowan.mitomotifs import seq2sites
from oldowan.mitomotifs import sites2str
from oldowan.mitomotifs import sites2seq

from oldowan.fasta import fasta

DLOOPplus = range(15800,16570) + range(1,1500)

dloops_fn   = os.path.join(here, 'dloops.fasta')
ofn         = os.path.join(here, 'fail_dloop.txt')
fail_fn     = os.path.join(here, 'dloop_expect_to_fail.txt')

expect_to_fail = []
for line in open(fail_fn, 'U'):
    expect_to_fail.append(int(line.strip()[:-1]))

of = open(ofn, 'w')

count = 0
for entry in fasta(dloops_fn):
    count += 1
    if count not in expect_to_fail:
        try:
            sites = seq2sites(entry["sequence"])
            seq = sites2seq(sites, region=DLOOPplus)
            seq = seq.replace('-', '')
            if entry["sequence"] in seq:
                print count, entry["name"], sites2str(sites)
            else:
                print 'FAILED', count, entry["name"]
                of.write('%d, %s\n' % (count, entry["name"]))
        except Exception, e:
            print 'FAILED',e, count, entry["name"]
            of.write('%d, %s\n' % (count, entry["name"]))
            of.flush()
of.close()
