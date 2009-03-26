"""This is the oldowan.mtconvert package."""

import os

VERSION = open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'VERSION')).read().strip()

__all__ = ['sites2seq', 'seq2sites', 'sites2str'] 

# Typically, this module will be installed as a sub-package under
# the oldowan package namespace. In this case, the primary import
# will work. 
#
# However, I have also included this sub-package directly in other
# applications. In this case, the namespace import approach fails and
# we fall back to direct import.
try:
    from oldowan.mtconvert.sites2seq import sites2seq
    from oldowan.mtconvert.str2sites import str2sites
    from oldowan.mtconvert.seq2sites import seq2sites
    from oldowan.mtconvert.sites2str import sites2str
except:
    from sites2seq import sites2seq
    from str2sites import str2sites
    from seq2sites import seq2sites
    from sites2str import sites2str
