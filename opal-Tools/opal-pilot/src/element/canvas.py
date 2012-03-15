from PyQt4.QtGui import QSizePolicy
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib import pyplot

class PlotCanvas(FigureCanvas):
    """
    Element responsible for rendering the plots
    """
    
    x = []
    y = []
    labels = {
        'x': 'undefined', 
        'y': 'undefined'
    }
    axes = None
    grid_state = True
    
    def __init__(self, parent=None, width=5, height=4, dpi=72):
        figure = pyplot.figure(dpi=dpi, facecolor='w', edgecolor='k')
        figure.subplots_adjust(left=0.18)
        
        self.axes = figure.add_subplot(111, adjustable='datalim')
        
        FigureCanvas.__init__(self, figure)
        self.setParent(parent)
        
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
               
    def grid(self, show=None):    
        if show is not None:
            self.grid_state = bool(show)
        self.axes.grid(self.grid_state)
        
    def setLabel(self, x=None, y=None):
        if x is not None:
            self.axes.set_xlabel(x)
        if y is not None:
            self.axes.set_ylabel(y)
        
    def addPlot(self, x, y, label):
        self.axes.plot(x, y, label=label)
        
    def legend(self):
        self.axes.legend(loc='best', prop={'size': 'x-small'}, fancybox=True)
        
    def clear(self):
        self.axes.clear()
        
    def plot(self, x=None, y=None, x_label=None, y_label=None):       
        if x is not None:
            self.x = x
        if y is not None:
            self.y = y
        if x_label is not None:
            self.labels['x'] = x_label
        if y_label is not None:
            self.labels['y'] = y_label

        self.axes.clear()
        for axis in ['x', 'y']:
            self.setLabel(axis, self.labels[axis])      
        self.axes.plot(self.x, self.y)
        self.axes.autoscale()
        self.grid()        
        pyplot.show()