# from matplotlib import pyplot
from PyQt4.QtGui import *
from PyQt4.QtCore import *

from . import View
from model.statistic import StatisticColumn
from element.canvas import PlotCanvas


class PlotView(View):
        
    def __init__(self):
        View.__init__(self)
    
    def setupUi(self):
        self.setWindowTitle('OPAL Pilot > Plotter')
        
        self.canvas = PlotCanvas()
        self.combo_x = QComboBox()
        self.combo_y = QComboBox()
        self.checkbox_grid = QCheckBox()
        self.checkbox_grid.setChecked(True)
        self.checkbox_legend = QCheckBox()
        self.checkbox_legend.setChecked(True)
        self.save_plot = QPushButton('Print', self)
        self.x_start = QLineEdit()
        self.x_end = QLineEdit()
        label_x = QLabel('x Axis')
        label_y = QLabel('y Axis')
        label_grid = QLabel('Show grid')
        label_legend = QLabel('Show legend')
        label_start = QLabel('Scale start')
        label_end = QLabel('Scale end')
        label_x.setBuddy(self.combo_x)
        label_y.setBuddy(self.combo_y)
        label_grid.setBuddy(self.checkbox_grid)
        label_legend.setBuddy(self.checkbox_legend)
        label_start.setBuddy(self.x_start)
        label_end.setBuddy(self.x_end)  
        
        layout_options = QGridLayout()
        layout_options.addWidget(label_start, 0, 0)
        layout_options.addWidget(self.x_start, 0, 1)
        layout_options.addWidget(label_end, 1, 0)
        layout_options.addWidget(self.x_end, 1, 1)
        layout_options.addWidget(label_x, 2, 0)
        layout_options.addWidget(self.combo_x, 2, 1)
        layout_options.addWidget(label_y, 3, 0)
        layout_options.addWidget(self.combo_y, 3, 1)
        layout_options.addWidget(label_grid, 4, 0)
        layout_options.addWidget(self.checkbox_grid, 4, 1)
        layout_options.addWidget(label_legend, 5, 0)
        layout_options.addWidget(self.checkbox_legend, 5, 1)
        layout_options.addWidget(self.save_plot, 6, 0)
        
        layout_main = QVBoxLayout()
        layout_main.addWidget(self.canvas)
        layout_main.addLayout(layout_options)       
        self.setLayout(layout_main)
        
    def addColumn(self, description, axis=StatisticColumn.AXIS_Y):
        """Add a column to the x or y combo box and return its index."""
        if axis == StatisticColumn.AXIS_X:
            combo = self.combo_x
        else:
            combo = self.combo_y
        count = combo.count()
        combo.addItem(description)
        return count


class InputMaskView(View):
    
    """
    Create view for the input mask.
    
    Each parameter group is rendered in its own tab. Each parameter consists
    of a label, a input box and optionally a unit.
    """
    
    def __init__(self):
        View.__init__(self)
        self.__tabbar = None
        self.__tabs = {}
        self.__parameters = {}      # Mapping of parameter_id to QLineEdit
        self.__servers = {}         # Mapping of QComboBox index to server_id
        self.__SERVER_GROUP_ID = 0
        self.__NOTE_GROUP_ID = -1
    
    def setupUi(self):
        """ Tab bar """
        self.__tabbar = QTabWidget()
        self.__tabbar.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)        
        self.__tabbar.setAutoFillBackground(False)
        self.__tabbar.setTabShape(QTabWidget.Rounded)
        self.__tabbar.setUsesScrollButtons(True)
        self.__tabbar.setDocumentMode(True)
        self.__tabbar.setTabsClosable(False)
        self.__tabbar.setMovable(False)
        
        """ Server Tab """
        self.addGroup(self.__NOTE_GROUP_ID, 'Note')
        self.addGroup(self.__SERVER_GROUP_ID, 'Server')
        self.__servercombo = QComboBox()
        self.__cpus = QLineEdit()
        self.__cpus.setText(str(8))
        self.__cpus.setValidator(QIntValidator())
        self.__note_title = QLineEdit()
        self.__note = QTextEdit()
        layout = self.__tabs[self.__NOTE_GROUP_ID]['layout']    # set in addGroup
        layout.addWidget(QLabel('Title'), 0, 0)
        layout.addWidget(self.__note_title, 0, 1)
        layout.addWidget(QLabel('Note'), 1, 0)
        layout.addWidget(self.__note, 1, 1)
        layout = self.__tabs[self.__SERVER_GROUP_ID]['layout']
        layout.addWidget(QLabel('Server'), 0, 0)
        layout.addWidget(self.__servercombo, 0, 1)
        layout.addWidget(QLabel('CPU Cores'), 1, 0)
        layout.addWidget(self.__cpus, 1, 1)       
        
        """ OK/Cancel buttons """
        self.buttons = QDialogButtonBox()
        self.buttons.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        
        """ Main layout """
        vbox = QVBoxLayout(self)
        vbox.setSpacing(0)
        vbox.setMargin(9)

        BoxScroll = QScrollArea()
        BoxScroll.setWidgetResizable(True)
        BoxScroll.setWidget(self.__tabbar)
        vbox.addWidget(BoxScroll)
        
        #vbox.addWidget(self.__tabbar)
        vbox.addWidget(self.buttons)
        
    def addGroup(self, id, alias):
        """ Create a new tab. """
        widget = QWidget()
        layout = QGridLayout(widget)
        self.__tabbar.addTab(widget, alias)
        self.__tabs[id] = {'layout': layout, 'line': 0}
        
    def addParameter(self, group_id, id, value, ui_label, ui_unit, ui_tooltip):
        """ Add parameter and return the QLineEdit instance """       
        layout = self.__tabs[group_id]['layout']
        line = self.__tabs[group_id]['line']
        
        if ui_unit is None:
            ui_unit = ''
        
        edit = QLineEdit()
        if value is not None:
            edit.setText(value)
        label = QLabel()
        label.setText(ui_label)
        label.setBuddy(edit)
        if ui_tooltip is not None:
            label.setToolTip(ui_tooltip)
        unit = QLabel()
        unit.setText(ui_unit)
        
        layout.addWidget(label, line, 0)
        layout.addWidget(edit, line, 1)
        layout.addWidget(unit, line, 2)
        
        self.__tabs[group_id] = {'layout': layout, 'line': line+1}
        return edit
        
    def addServer(self, id, ui_alias):
        combo_index = self.__servercombo.count()
        self.__servers[combo_index] = id
        self.__servercombo.addItem(ui_alias)
        return combo_index
    
    def server(self):
        return self.__servercombo.currentIndex()
    
    def serverCores(self):
        return self.__cpus.text()
    
    def setDefault(self, cores=None, note_title=None, note=None):
        if cores is not None:
            self.__cpus.setText(str(cores))
        if note_title is not None:
            self.__note_title.setText(str(note_title))
        if note is not None:
            self.__note.setText(str(note))
    
    def note(self):
        note_title = self.__note_title.text()
        note = self.__note.toPlainText()
        print note
        return (note_title, note)
        
        
class VisibleColumnPickerView(View):
    
    def __init__(self):
        View.__init__(self)
        self.__groups = {}
        
        self.__line = 0
        self.__column = 0
        self.__columns = 2
        
    def setupUi(self):
        self.__layout = QVBoxLayout(self)
        self.__options = QGridLayout()
        
        """ OK/Cancel buttons """
        self.buttons = QDialogButtonBox()
        self.buttons.setStandardButtons(QDialogButtonBox.Ok)
        
        self.__layout.addLayout(self.__options)
        self.__layout.addWidget(self.buttons)                
    
    def addGroup(self, group_id, ui_label, checked=True):
        layout = QVBoxLayout()
        group = QGroupBox()
        group.setCheckable(True)
        group.setChecked(checked)
        group.setFlat(True)
        group.setTitle(ui_label)
        group.setLayout(layout)
        self.__groups[group_id] = layout
        self.__options.addWidget(group, self.__line, self.__column)
        
        """ x-Column Layout """
        if self.__column < self.__columns-1:
            self.__column += 1
        else:
            self.__line += 1
            self.__column = 0
        
        return group
        
    def addColumn(self, group_id, ui_label, checked=True):
        checkbox = QCheckBox()
        checkbox.setCheckState(checked)
        checkbox.setTristate(False)
        checkbox.setText(ui_label)
        self.__groups[group_id].addWidget(checkbox)
        
        return checkbox
