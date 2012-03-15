from time import time
from functools import partial

from PyQt4.QtCore import QSize, QTimer, SIGNAL
from PyQt4.QtGui import QWidget, QVBoxLayout, QTabWidget, QSizePolicy, QSpacerItem, QColor, QHBoxLayout, QTreeWidget, QTreeWidgetItem, QAbstractItemView, QCheckBox, QMenuBar

from . import View
from element.button import AcceleratorButton

class MainWindowView(View):
    
    """
    Create a view for the main application window.
    
    The view renders each accelerator in its own tab.
    """
    
    def __init__(self):
        View.__init__(self)        
        self.central_widget = None
        self.__central_layout = None
        self.__tabbar = None
    
    def setupUi(self):
        self.setWindowTitle('OPAL Pilot')
        
        # Base frame      
        self.central_widget = QWidget()
        self.__central_layout = QVBoxLayout(self.central_widget)
        self.__central_layout.setSpacing(0)
        self.__central_layout.setMargin(0)
        
        # Menu bar
        menubar = QMenuBar(self.central_widget)
        file_menu = menubar.addMenu('File')
        file_menu.addAction('Exit', self.close)  
        
        # Tab bar holding accelerators
        self.tabbar = QTabWidget(self.central_widget)
        self.tabbar.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)        
        self.tabbar.setAutoFillBackground(False)
        self.tabbar.setTabShape(QTabWidget.Rounded)
        self.tabbar.setUsesScrollButtons(True)
        self.tabbar.setDocumentMode(True)
        self.tabbar.setTabsClosable(False)
        self.tabbar.setMovable(False)
        
        self.__central_layout.addWidget(menubar)
        self.__central_layout.addWidget(self.tabbar)
        
    def addAccelerator(self, accelerator_controller, ui_alias):
        self.tabbar.addTab(accelerator_controller.widget, ui_alias)     


class AcceleratorView(View):
    
    #    __HIGHLIGHT_COLOR_ON = QColor(244, 244, 244) # PSI light-grey
    __HIGHLIGHT_COLOR_ON = QColor(184, 215, 239) 
    __HIGHLIGHT_COLOR_OFF = QColor(255, 255, 255) # White
    
    COLUMNID_STATE = -4    
    
    def __init__(self):
        View.__init__(self)
        self.widget = None
        self.__button_box = None
        self.__checkbox_box = None
        self.__listcheckbox = []
        self.__spacer = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
        self.__button_size = QSize(50, 40)
        self.__column_order = []
        self.__highlighted = [] # [(timestamp_start, QTreeWidgetItem), ...]
        self.addTimer(10000, self.__timerTurnOffHighlights)
    
    def setupUi(self):
        self.widget = QWidget()
        vert_layout = QVBoxLayout(self.widget)
        self.__button_box = QHBoxLayout()
        self.__button_box.addItem(self.__spacer)       
        self.tree = QTreeWidget()
        self.tree.setDragEnabled(True)
        self.tree.setItemsExpandable(False)
        self.tree.setRootIsDecorated(False)
        self.tree.setSortingEnabled(True)
        # self.tree.setDragDropMode(QAbstractItemView.)
        self.tree.setSelectionMode(QAbstractItemView.ExtendedSelection)
        vert_layout.addLayout(self.__button_box)
        vert_layout.addWidget(self.tree)
        
    def addButton(self, ui_label, drop_event_type=None, drop_policy=AcceleratorButton.DROPS_ACCEPT_NONE):
        """ Add a button to the button box and return the widget """
        button = AcceleratorButton(self, drop_event_type, drop_policy)
        button.setText(ui_label)
        button.setMinimumSize(self.__button_size)        
        self.__button_box.removeItem(self.__spacer)
        self.__button_box.addWidget(button)
        self.__button_box.addItem(self.__spacer)
        return button
    
    def addColumn(self, parameter_id, ui_label, ui_tooltip=None, first=False, hidden=False):
        self.__column_order.append(parameter_id)
        header_item = self.tree.headerItem()
        # The header item is always created with one column, so automatic 
        # detection would fail in that case.
        if first is False:       
            index = header_item.columnCount()
        else:
            index = 0
        header_item.setText(index, ui_label)
        # self.addCheckBox(index, ui_label)
        if ui_tooltip is not None:
            header_item.setToolTip(index, ui_tooltip)            
        return index
            
    def addRow(self, contents, highlight=False):
        """ Add row to the TreeWidget and return TreeWidgetItem. """
        elements = []
        for parameter_id in self.__column_order:
            if contents.has_key(parameter_id):
                elements.append(contents[parameter_id])
            else:
                elements.append('')
        item = QTreeWidgetItem(elements)
        self.tree.addTopLevelItem(item)
        
        if highlight is True:
            self.__highlight(item)
       
        return item
    
    def editRow(self, item, column_id, new_value):
        column = self.__column_order.index(column_id)
        item.setText(column, new_value)
        self.__highlight(item)
    
    def __highlight(self, widget_item):
#        for i in range(widget_item.columnCount()):
#            widget_item.setBackgroundColor(i, self.__HIGHLIGHT_COLOR_ON)
        self.__colorTreeWidgetItem(widget_item, self.__HIGHLIGHT_COLOR_ON)
        self.__highlighted.append((int(time()), widget_item))
            
    def __timerTurnOffHighlights(self):
        for tuple in self.__highlighted:
            time_started, item = tuple
            if time_started < int(time())-20:
                self.__colorTreeWidgetItem(item, self.__HIGHLIGHT_COLOR_OFF)
                self.__highlighted.remove(tuple)
                
    def __colorTreeWidgetItem(self, item, color):
        for i in range(item.columnCount()):
            item.setBackgroundColor(i, color)
        
