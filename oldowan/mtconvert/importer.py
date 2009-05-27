from oldowan.mtconvert.popset import Sample, Population, PopSet
from oldowan.mtconvert.coverage import Coverage
from oldowan.mtconvert.str2sites import str2sites
from oldowan.mtconvert.error import MtconvertError

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
             doi            = False,
             pmid           = False,
             ):
    """Load mitochondrial haplotype definitions from csv file.
    
    Column numbers are to be provided in python-standard 0-based counting.
    """

    if population is False:
        population = "Unknown"

    errors = []
    line_number = header

    if sites:
        if len(sites) != len(sites_on_rCRS):
            errors.append((0, MtconvertError("When sites are included, sites_on_rCRS must match")))

    # base coverage is given by hvr1 and hvr2 covers
    segments = list(x for x in [hvr1_covers, hvr2_covers] if x)

    # every sample needs an id,
    #  if it is given in the file, we use that
    #  otherwise create one in the format 's#'
    if sample_id is not False:
        def sample_generator(line):
            sids = line[sample_id].split(sample_id_sep)
            for sid in sids:
                yield (sid, population)
    else:
        def sample_generator(line):
            global s
            count = 0
            pop_names = []
            if pop_with_n:
                max = sum(num(line[x]) for x in n)
                for i in range(len(n)):
                    pop_names = pop_names + [population[i]] * num(line[n[i]])
            elif n:
                max = num(line[n])
                pop_names = [population] * max
            else:
                max = 1
                pop_names = [population]
            while count < max:
                pop = pop_names[count]
                count += 1
                s += 1
                yield ("s%d" % s, pop)

    # start the reader
    reader = csv.reader(open(file, 'rU'))
    samples = []
    for l in list(x for x in reader)[header:]:

        line_number += 1

        ###################################################
        # Read in polymorphism data
        ###################################################

        polymorphisms = []

        this_sample_segments = list(x for x in segments)

        if hvr1:
            try:
                polys = str2sites(l[hvr1], add16k=add16k)
                polymorphisms += polys
            except MtconvertError, e:
                errors.append( (line_number, e) )

        if hvr2:
            try:
                polys = str2sites(l[hvr2])
                polymorphisms += polys
            except MtconvertError, e:
                errors.append( (line_number, e) )

        if sites:
            for i in range(0,len(sites)):
                site_index = sites[i]
                position   = sites_on_rCRS[i]
                value      = l[site_index].strip().upper()
                if value in ('A','G','C','T'):
                    try:
                        poly = str2sites('%d%s' % (position,value))
                        polymorphisms.append(poly)
                        this_sample_segments.append(position)
                    except MtconvertError, e:
                        errors.append( (line_number, e) )

        coverage = Coverage(*this_sample_segments)

        ###################################################
        # Extract haplogroup
        ###################################################

        hap = None
        if haplogroup is not False:
            hap = l[haplogroup].strip()

        sample_ids = sample_generator(l)        
        for (sid,pop) in sample_ids:
            samples.append(Sample(id=sid,
                                  haplogroup=hap,
                                  population=pop,
                                  coverage=coverage,
                                  polymorphisms=polymorphisms))

    pop = Population(samples=samples)
    return PopSet(populations=[pop],
                  errors=errors)

