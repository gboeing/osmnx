
class UnknownNetworkType(ValueError):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)


class InsufficientNetworkQueryArguments(ValueError):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)


class EmptyOverpassResponse(ValueError):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)


class InvalidDistanceType(ValueError):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)


