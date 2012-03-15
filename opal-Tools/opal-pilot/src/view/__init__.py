import uuid
import warnings
from functools import partial

from PyQt4.QtCore import QTimer, SIGNAL

class View(object):
    
    def __init__(self):
        self.__timers = {}
    
    def setupUi(self):
        """Set up the UI elements.""" 
        pass
        
    def addTimer(self, interval, callable, first_run=None, label=None):
        """Add a timer that runs in fixed intervals.
        
        The timer will be first issued after the first interval has passed.
        Timings: (interval, interval*2, interval*3, ...).
        
        Optionally the timer can be started after first_run milliseconds for
        the first time. The intervals will then start after the first run.
        Timings: (first_run, first_run+interval, first_run+interval*2, ...)
        
        Keyword arguments:
        interval -- the interval in milliseconds
        callable -- a python callable
        first_run -- the time after that the callable is issued for the 
                     first time.
         
        """
        if label is None:
            label = uuid.uuid4()
        if first_run is not None:
            QTimer.singleShot(first_run, partial(self.addTimer, interval=interval, callable=callable, label=label))
            QTimer.singleShot(first_run, callable)
        else:
            timer = QTimer()
            timer.start(interval)
            self.connect(timer, SIGNAL('timeout()'), callable)
            self.__timers[label] = timer
            
    def stopTimer(self, label):
        """Stop a timer."""
        if label in self.__timers.keys():
            self.__timers[label].stop()
            del self.__timers[label]
        else:
            warnings.warn('There is no timer with the label: %s' % label, Warning, stacklevel=2)