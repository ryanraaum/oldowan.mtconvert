from oldowan.mtconvert.popset import Sample, Population, PopSet
from oldowan.mtconvert.coverage import Coverage

import csv

# global counter variable for newly-created sample ids
s = 0

def num(x):
    if x == '':
        return 0
    else:
        return int(x)

def load_csv(file,
             header         = 1,     # number of rows to skip for header info
             hvr1           = False, # column number if present
             hvr1_covers    = [16024,16365],
             add16k         = True,  # add 16000 to every hvr1 site?
             hvr2           = False, # column number if present
             hvr2_covers    = [73,340],
             rflps          = False, # column number if present
             rflp_format    = False, # ? how to implement this ?
             sites          = False, # column number(s) if present
             sites_on_rCRS  = [], # matched entry or list for sites columns
             haplogroup     = False, # column number if present
             sample_id      = False, # column number if present
             sample_id_sep  = ',',   # what separates multiple ids?
             haplotype_id   = False, # column number if present
             pop_with_n     = False, # are N's arranged by population?
             n              = False, # column number(s) if present
             population     = False, # column number or name or if pop_with_n is True, 
                                     # names to go with N columns 
             ):
    """Load mitochondrial haplotype definitions from csv file."""

    # calculate coverage
    segments = list(x for x in [hvr1_covers, hvr2_covers] if x)
    if sites_on_rCRS:
        segments = segments + sites_on_rCRS
    coverage = Coverage(*segments)

    # every sample needs an id,
    #  if it is given in the file, we use that
    #  otherwise create one in the format 's#'
    if sample_id:
        def sample_generator(line):
            return line[sample_id-1].split(sample_id_sep)
    else:
        def sample_generator(line):
            global s
            count = 0
            if pop_with_n:
                max = sum(num(line[x-1]) for x in n)
            elif n:
                max = num(line[n-1])
            else:
                max = 1
            while count < max:
                count += 1
                s += 1
                yield ["s%d" % s]

    # start the reader
    reader = csv.reader(open(file, 'rU'))
    samples = []
    for l in list(x for x in reader)[header:]:
        sample_ids = sample_generator(l)        
        for sid in sample_ids:
            samples.append(Sample(id=sid))

    pop = Population(coverage=coverage, samples=samples)
    return PopSet(populations=[pop])
