import warnings

def deprecated(func):
    """
    Decorator to mark functions deprecated.
    
    The Decorator can be used to mark functions deprecated. A warning will be
    issued when the function is used.
    
    Source: http://code.activestate.com/recipes/391367-deprecated/
    """
    def newFunc(*args, **kwargs):
        warnings.warn("Call to deprecated function %s." % func.__name__,
                      category=DeprecationWarning, stacklevel=2)
        return func(*args, **kwargs)
    newFunc.__name__ = func.__name__
    newFunc.__doc__ = func.__doc__
    newFunc.__dict__.update(func.__dict__)
    return newFunc

hotshotProfilers = {}
def hotshotit(func):
    """
    Profile pure-Python code, using hotshot.
     
    Source: http://code.activestate.com/recipes/576656-quick-python-profiling-with-hotshot/    
    """
    def wrapper(*args, **kw):
        import hotshot
        global hotshotProfilers
        prof_name = func.func_name+".prof"
        profiler = hotshotProfilers.get(prof_name)
        if profiler is None:
            profiler = hotshot.Profile(prof_name)
            hotshotProfilers[prof_name] = profiler
        return profiler.runcall(func, *args, **kw)
    return wrapper