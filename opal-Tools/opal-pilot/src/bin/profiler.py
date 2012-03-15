"""
Display results from hotshot profiler files (.prof)

Usage:
python profiler.py <filename>.prof

Example of how to create a profile:
    ...
    from util.decorator import hotshotit
    
    class SomeClass():
    
        @hotshotit                <--
        def someMethod(self):
    ...
"""
import sys
import hotshot, hotshot.stats
stats = hotshot.stats.load(sys.argv[1])
stats.strip_dirs()
stats.sort_stats('time', 'calls')
stats.print_stats(25)
