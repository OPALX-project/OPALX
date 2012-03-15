class Singleton(type):
    """
    Singleton pattern taken from Wikipedia
    http://en.wikipedia.org/wiki/Singleton_pattern#Python
    
    Usage:
    class MyClass(object):
        __metaclass__ = Singleton
    
    """
    
    def __init__(cls, name, bases, dict):
        super(Singleton, cls).__init__(name, bases, dict)
        cls.instance = None
 
    def __call__(cls, *args, **kw):
        if cls.instance is None:
            cls.instance = super(Singleton, cls).__call__(*args, **kw)
 
        return cls.instance