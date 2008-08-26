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

dloops_fn = os.path.join(here, 'dloops.fasta')
ofn = os.path.join(here, 'fail_dloop.txt')
of = open(ofn, 'w')

count = 0
for entry in fasta(dloops_fn):
    count += 1
    if count > 20:
	    break
    try:
        sites = seq2sites(entry["sequence"])
        seq = sites2seq(sites, what=DLOOPplus)
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