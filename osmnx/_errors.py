################################################################################
# Module: errors.py
# Description: Custom errors
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

class EmptyOverpassResponse(ValueError): # pragma: no cover
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)


class InsufficientNetworkQueryArguments(ValueError): # pragma: no cover
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)


class InvalidDistanceType(ValueError): # pragma: no cover
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)


class UnknownNetworkType(ValueError): # pragma: no cover
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)
