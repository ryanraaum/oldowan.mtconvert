
class MtconvertError(Exception):
    """Exception raised for errors in the mtconvert module.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message

