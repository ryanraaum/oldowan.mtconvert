from oldowan.mtconvert.coverage import Coverage

class Sample(object):

    def __init__(self, id=None,
                       haplotype_id=None,
                       population_id=None,
                       coverage=Coverage((16024,16365)),
                       doi=None,
                       pmid=None,
                       ):
        self._id                = id
        self._haplotype_id      = haplotype_id
        self._population_id     = population_id
        self._coverage          = coverage
        self._doi               = doi
        self._pmid              = pmid


class Population(object):

    def __init__(self, samples=[],
                       coverage=Coverage((16024,16365)),
                       ):

        self.__samples = samples
        self.__coverage = coverage

    def __get_num_samples(self):
        return len(self.__samples)

    num_samples = property(fget=__get_num_samples)

    def __get_coverage(self):
        return self.__coverage

    coverage = property(fget=__get_coverage)


class PopSet(object):

    def __init__(self, populations=[]):
        self.__populations = populations
        if not populations:
            self.__coverage = Coverage()
        else:
            c = populations[0].coverage
            def collect(x):
                c = c.intersection(x)
            [collect(x) for x in populations[1:]]
            self.__coverage = c

    def __get_num_populations(self):
        return len(self.__populations)

    num_populations = property(fget=__get_num_populations)

    def __get_num_samples(self):
        return sum(list(x.num_samples for x in self.__populations)) 

    num_samples = property(fget=__get_num_samples)

    def __get_coverage(self):
        return self.__coverage

    coverage = property(fget=__get_coverage)
