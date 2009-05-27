from oldowan.mtconvert.coverage import Coverage


class Sample(object):

    def __init__(self, id,
                       haplotype_id=None,
                       haplogroup=None,
                       population=None,
                       coverage=Coverage(),
                       polymorphisms=[],
                       doi=None,
                       pmid=None,
                       ):
        self.__id                = id
        self.__haplotype_id      = haplotype_id
        self.__haplogroup        = haplogroup
        self.__population        = population
        self.__coverage          = coverage
        self.__polymorphisms     = polymorphisms
        self.__doi               = doi
        self.__pmid              = pmid

    def __get_id(self):
        return self.__id

    id = property(fget=__get_id)

    def __get_haplogroup(self):
        return self.__haplogroup

    haplogroup = property(fget=__get_haplogroup)

    def __get_coverage(self):
        return self.__coverage

    coverage = property(fget=__get_coverage)

    def __get_population(self):
        return self.__population

    population = property(fget=__get_population)


class Population(object):

    def __init__(self, samples=[],
                       ):

        self.__samples = samples
        if not samples:
            self.__coverage = Coverage()
        else:
            c = samples[0].coverage
            for x in samples[1:]:
                c = c.intersection(x.coverage)
            self.__coverage = c

        self.__samples_dict = {}
        for s in self.__samples:
            # TODO: this will break if more than one sample has the same id
            self.__samples_dict[s.id] = s


    def __get_num_samples(self):
        return len(self.__samples)

    num_samples = property(fget=__get_num_samples)

    def __get_coverage(self):
        return self.__coverage

    coverage = property(fget=__get_coverage)

    def __get_samples(self):
        return self.__samples

    samples = property(fget=__get_samples)

    def sample_by_id(self, id):
        return self.__samples_dict[id]

class PopSet(object):

    def __init__(self, populations=[], errors=None):
        self.__populations = populations
        if not populations:
            self.__coverage = Coverage()
        else:
            c = populations[0].coverage
            for x in populations[1:]:
                c = c.intersection(x.coverage)
            self.__coverage = c
        self.__errors = errors

    def __get_num_populations(self):
        return len(self.__populations)

    num_populations = property(fget=__get_num_populations)

    def __get_populations(self):
        return self.__populations

    populations = property(fget=__get_populations)

    def __get_num_samples(self):
        return sum(list(x.num_samples for x in self.__populations)) 

    num_samples = property(fget=__get_num_samples)

    def __get_coverage(self):
        return self.__coverage

    coverage = property(fget=__get_coverage)

    def __get_errors(self):
        return self.__errors

    errors = property(fget=__get_errors)

