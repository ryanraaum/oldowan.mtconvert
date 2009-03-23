from types import StringType, ListType
from oldowan.polymorphism import Polymorphism

def sites2str(sites):
    """Transform a list of Polymorphisms to a string.

    """
    processing = []
    for x in sites:
        # unnested (unambiguous) sites get directly pushed to the string
        if isinstance(x,Polymorphism):
            processing.append(str(x))
        # ambiguous sites are nested in lists
        elif type(x) == ListType:
            current = []
            for y in x:
                if type(y) == ListType:
                    if len(y) == 1:
                        current.append(str(y[0]))
                    elif len(y) > 1:
                        as_str = list(str(x) for x in y)
                        interior = ' '.join(as_str)
                        current.append('(' + interior + ')')
                    else:
                        raise Exception('format error in sites2str')
                else:
                    raise Exception('format error in sites2str')
            interior = ' or '.join(current)
            processing.append('(' + interior + ')')
        else:
            raise Exception('format error in sites2str')
    return ' '.join(processing)

