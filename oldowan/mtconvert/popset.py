from oldowan.mtconvert.coverage import Coverage

class Population(object):

    def __init__(self):
        self.__samples = []
        self.__coverage = Coverage((16024,16365))

    def __get_num_samples(self):
        return len(self.__samples)

    num_samples = property(fget=__get_num_samples)

    def __get_coverage(self):
        return self.__coverage

    coverage = property(fget=__get_coverage)


class PopSet(object):

    def __init__(self):
        self.__populations = []
        self.__coverage = Coverage((16024,16365))

    def __get_num_populations(self):
        return len(self.__populations)

    num_populations = property(fget=__get_num_populations)

    def __get_num_samples(self):
        return sum(list(x.num_samples for x in self.__populations)) 

    num_samples = property(fget=__get_num_samples)

    def __get_coverage(self):
        return self.__coverage

    coverage = property(fget=__get_coverage)
