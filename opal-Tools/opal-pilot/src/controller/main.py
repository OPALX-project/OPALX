import os
import time
from functools import partial

from PyQt4 import Qt
from PyQt4.QtGui import QMainWindow, QDesktopWidget, QWidget
from PyQt4.QtCore import SIGNAL, QSettings, QSize, QPoint

from model.accelerator import Accelerator
from model.simulation import Simulation
from model.statistic import Statistic
from view.main import MainWindowView, AcceleratorView
from controller.dialog import InputMaskController, PlotController, VisibleColumnPickerController
from element.button import AcceleratorButton
from element.thread import Postprocessing, PollServer
from element.canvas import PlotCanvas       

class MainWindowController(QMainWindow, MainWindowView):
    """Render main window.
    
    @todo Add a menu bar with exit button
    @todo Add a status bar
    """
    
    def __init__(self):
        QMainWindow.__init__(self)
        MainWindowView.__init__(self)
        
        self.__accelerators = {}        # {accelerator_id: AcceleratorController, ...}
        self.__threads = [Postprocessing(), PollServer()]  
        accelerator = Accelerator()
        
        self.setupUi()        
        
        for acc in accelerator.names():
            controller = AcceleratorController(acc['id'])
            self.addAccelerator(controller, acc['ui_alias'])
            self.__accelerators[acc['id']] = controller
            
        desktop = QDesktopWidget()
        screen = desktop.screenGeometry()
                             
        self.setCentralWidget(self.central_widget)
        self.readSettings()
        self.show()
        
        self.connect(self.tabbar, SIGNAL('currentChanged(int)'), self.tabChanged)
        self.addTimer(label='threads', interval=10000, first_run=5000, callable=self.__timerThreads)
        
    def __timerThreads(self):
        threads = self.__threads
        for i in range(len(threads)):
            threads[i]
            if not threads[i].is_alive():
                cls = threads[i].__class__
                threads[i] = cls()
                threads[i].start()
                
    def writeSettings(self):
        settings = QSettings()
        settings.beginGroup('MainWindow')
        settings.setValue('size', self.size())
        settings.setValue('position', self.pos())
        settings.endGroup()
        
    def readSettings(self):
        settings = QSettings()
        
        desktop = QDesktopWidget()
        screen = desktop.screenGeometry()
        
        settings.beginGroup('MainWindow')
        size = settings.value('size', QSize(screen.width(), screen.height()/2)).toSize()
        self.resize(size)
        geometry = self.geometry()
        position = settings.value('position', QPoint((screen.width()-geometry.width())/2, 0)).toPoint()
        self.move(position)
        accelerator_index = settings.value('accelerator_active', 0).toInt() # no idea why this returns a tuple
        self.tabbar.setCurrentIndex(accelerator_index[0])
        settings.endGroup()
        
    def closeEvent(self, event):
        # Save settings
        # MainWindow
        self.writeSettings()
        # Accelerators
        for accelerator in self.__accelerators.itervalues():
            accelerator.closeEvent()
        # Wait for threads to finish
        print "waiting for background threads to finish..."
        self.stopTimer('threads')
        while True:
            threads = self.__threads
            finished = True
            for thread in threads:
                if thread.is_alive():
                    finished = False
            if finished:
                break
            else:
                time.sleep(1)
        """ TODO: clean up tmp/ dir"""
        # self.__cleanTmp()
                                
        event.accept()
        
    def tabChanged(self, index):
        settings = QSettings()
        settings.setValue('MainWindow/accelerator_active', index)
        
    def __cleanTmp(self, path='tmp/', recursive=True, remove_folder=False):
        contents = os.listdir(path)
        for item in contents:
            if os.path.isfile(path + item):
                os.remove(path + item)
            elif recursive and item != '.svn' and os.path.isdir(path + item):
                self.__cleanTmp(path + item + "/", True, True)
        if remove_folder:
             os.rmdir(path)
             pass


class AcceleratorController(QWidget, AcceleratorView):
    
    """ TODO: Refactor the part of the code responsible for loading and applying the column visibility settings. """
    
    TIMER_UPDATE_SIMULATION = 10000 
    
    def __init__(self, accelerator_id):
        QWidget.__init__(self)
        AcceleratorView.__init__(self)
        
        self.__accelerator_id = accelerator_id
        self.__items = {}                        # {simulation_id: QTreeWidgetItem, ...}
        self.__columns = {}                      # {parameter_id: index, ...}
        self.__group_visibility = {}
        self.__parameter_visibility = {}        
        accelerator = Accelerator()
        simulation = Simulation()
        
        widget = self.setupUi()
        
        buttons = {
            1: ('New Simulation', AcceleratorButton.EVENT_NEW, AcceleratorButton.DROPS_ACCEPT),  
            2: ('Plot Results', AcceleratorButton.EVENT_PLOT, AcceleratorButton.DROPS_ACCEPT_MULTIPLE),
            3: ('Columns', None, AcceleratorButton.DROPS_ACCEPT_NONE), 
#           4: ('H5Root', None, AcceleratorButton.DROPS_ACCEPT_NONE),
#           5: ('Filter', None, AcceleratorButton.DROPS_ACCEPT_NONE),
        }
        for key in buttons.iterkeys():
            label, drop_event_type, drops_policy = buttons[key]
            button = self.addButton(label, drop_event_type, drops_policy)
            self.connect(button, SIGNAL('clicked(bool)'), partial(self.slotButtonClicked, key))
            
        # Visibility settings
        settings = QSettings()
        settings.beginGroup('AcceleratorVisibleColumns')
        size = settings.beginReadArray(self.__accelerator_id + '_groups')
        for i in range(size):
            settings.setArrayIndex(i)
            group_id = settings.value('group_id').toString()
            state = settings.value('state').toBool()
            self.__group_visibility[int(group_id)] = state
        settings.endArray()
        size = settings.beginReadArray(self.__accelerator_id + '_columns')
        for i in range(size):
            settings.setArrayIndex(i)
            parameter_id = settings.value('parameter_id').toString()
            state = settings.value('state').toBool()
            self.__parameter_visibility[int(parameter_id)] = state        
        settings.endArray()        
        settings.endGroup()     
            
        # Columns that are valid for all accelerators
        general_columns = [
            (-7, 'Date', True),
            (-6, 'Username', False),
            (-5, 'Note', False),
            (self.COLUMNID_STATE, 'State', False),
            (-3, 'ID', False),
            (-2, 'Process', False),
            (-1, 'Server', False),
            (0, 'Cores', False)
        ]
        for tuple in general_columns:
            id, label, first = tuple
            index = self.addColumn(id, label, first=first)
            self.__columns[id] = index
            hidden = self.isColumnHidden(0, id)
            self.toggleColumn(id, not hidden)
            
        # Add the columns that are specific for this accelerator.
        parameters = accelerator.parameters(accelerator_id)
        if parameters is not None:
            for parameter in parameters:
                parameter_id = int(parameter['id'])
                group_id = int(parameter['group_id'])
                index = self.addColumn(parameter_id, parameter['label'], parameter['tooltip'])
                self.__columns[parameter_id] = index
                hidden = self.isColumnHidden(group_id, parameter_id)
                self.toggleColumn(parameter_id, not hidden)
        
        # Fill the accelerator with data row
        simulation_ids = simulation.ids(accelerator_id)
        if len(simulation_ids) > 0:        
            simulations = simulation.values(simulation_ids)
            for id in simulation_ids:
                # state = simulations[id]['state']
                row = self.__prepareValuesForView(simulations[id], id)
                item = self.addRow(row)
                self.__items[id] = item
                      
        self.addTimer(self.TIMER_UPDATE_SIMULATION, self.__timerUpdateSimulations)
        
        self.readSettings()
        self.resizeColumnsToContents()
                        
        return widget
    
    def isColumnHidden(self, group_id, parameter_id):
        if group_id in self.__group_visibility.keys() and self.__group_visibility[group_id] == False:
            return True
        if parameter_id in self.__parameter_visibility.keys() and self.__parameter_visibility[parameter_id] == False:
            return True
        return False
        
    
    def __timerUpdateSimulations(self):
        simulation = Simulation()
        ids = simulation.ids(self.__accelerator_id, self.TIMER_UPDATE_SIMULATION)
        for id in ids:
            if id not in self.__items.keys():
                self.newRow(id)
            else:
                self.updateState(id)
                                    
    def __prepareValuesForView(self, simulation_values_item, simulation_id):
        x = simulation_values_item
        output = {
            -7: str(x['date']),
            -6: str(x['username']),
            -5: str(x['note_title']),
            -4: str(x['state']),
            -3: str(simulation_id),
            -2: str(x['process_id']),
            -1: str(x['server']),
             0: str(x['server_cores'])
        }
        for id in x['values'].iterkeys():
            output[int(id)] = x['values'][id] 
        return output
    
    def newRow(self, simulation_id):
        simulation = Simulation()
        list = simulation.values([simulation_id])
        if len(list) > 0:
            item = list[str(simulation_id)]
            row = self.__prepareValuesForView(item, simulation_id)
            widget = self.addRow(row, highlight=True)
            self.__items[str(simulation_id)] = widget
            
    def updateState(self, simulation_id):
        simulation = Simulation()
        state = simulation.getState(simulation_id)
        item = self.__items[simulation_id]
        self.editRow(item, self.COLUMNID_STATE, state)        
            
    def openInputMask(self, default_value_item=None):
        if default_value_item is None:
            mask = InputMaskController(self.__accelerator_id, simulation_id=None, parent=self)
        elif default_value_item in self.__items.values():
            for simulation_id in self.__items.keys():
                if self.__items[simulation_id] == default_value_item:
                    break
            mask = InputMaskController(self.__accelerator_id, simulation_id=simulation_id, parent=self)
        else:
            raise Exception('No such default_value_item')
        
    def openColumnPicker(self):
        picker = VisibleColumnPickerController(self.__accelerator_id, self)
        
    def resizeColumnsToContents(self):
        for i in range(self.tree.columnCount()):
            self.tree.resizeColumnToContents(i)
        self.tree.update()
        
    def toggleColumn(self, parameter_id, column_visible=None):
        index = self.__columns[parameter_id]
        if column_visible is not None:
            self.tree.setColumnHidden(index, not column_visible)            
        else:
            if self.tree.isColumnHidden(index):
                hide = False
            else:
                hide = True
            self.tree.setColumnHidden(index, hide)

        self.resizeColumnsToContents()
        
    def openPlottingWindow(self, items):
        simulation = Simulation()
        ids_dropped = []
        ids_finished = []
        for id in self.__items.iterkeys():
            if self.__items[id] in items:
                ids_dropped.append(id)
        for id in ids_dropped:
            state = simulation.getState(id)
            if state == simulation.STATE_FINISHED:
                ids_finished.append(id)
        if len(ids_finished) > 0:
            window = PlotController(parent=self, accelerator_id=self.__accelerator_id, simulation_ids=ids_finished)
                
    def slotButtonClicked(self, index):
        if index is 1:
            self.openInputMask()
        if index is 3:
            self.openColumnPicker()
            
    def readSettings(self):
        settings = QSettings()
        settings.beginGroup('Accelerator/'+self.__accelerator_id)
        sort_index, _ = settings.value('sort_column', 0).toInt()
        sort_order, _ = settings.value('sort_order', 0).toInt()
        self.tree.sortItems(sort_index, sort_order)
        settings.endGroup()
        
    def closeEvent(self):
        pass       
        
