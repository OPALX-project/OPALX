import sys
from functools import partial

from matplotlib import pyplot
from PyQt4.QtGui import QDialog
from PyQt4.QtCore import SIGNAL, SLOT, QSettings

from model.accelerator import Accelerator
from model.simulation import Simulation, SimulationStatisticColumn
from model.server import Server
from model.statistic import StatisticColumn
from view.dialog import InputMaskView, PlotView, VisibleColumnPickerView

from util.decorator import hotshotit

class PlotController(QDialog, PlotView):
    
    def __init__(self, accelerator_id, simulation_ids, parent=None):
        QDialog.__init__(self, parent)
        PlotView.__init__(self)
        
        self.__columns = {}
        self.__simulations = simulation_ids
        self.__simulation_statistic_column = SimulationStatisticColumn() # make use of the caching
        self.__index2id = {} # {combo_box_index: {axis: statistic_column_id}, ...}
        self.__current_plot = {
            StatisticColumn.AXIS_X: None, 
            StatisticColumn.AXIS_Y: None
        } # statistic_column_ids of the currently plotted axis
        self.__grid_state = True
        self.__legend_state = True
        self.__x_start = None
        self.__x_end = None
        
        accelerator = Accelerator()
        statistic_column = StatisticColumn()
        column_ids = self.__simulation_statistic_column.columnIds(simulation_ids)
        column_ids = [x[0] for x in column_ids]
        columns = statistic_column.columns(column_ids)
        for column in columns:
            id = column['id']
            del column['id']
            self.__columns[int(id)] = column
             
        self.setupUi()
        self.setWindowTitle('OPAL Pilot > %s > Plotter' % accelerator.alias(accelerator_id))
        
        keys = self.__columns.keys()
        keys.sort()
        for id in keys:
            column = self.__columns[id]
            # set the first row of every axis as default
            if self.__current_plot[column['ui_axis']] is None:
                self.__current_plot[column['ui_axis']] = id
            self.addColumn(id, column['description'], column['ui_axis'])
        self.plot() # plot with default axis
        self.show()
        
        self.connect(self.combo_x, SIGNAL('currentIndexChanged(int)'), partial(self.selectionChanged, axis=StatisticColumn.AXIS_X))
        self.connect(self.combo_y, SIGNAL('currentIndexChanged(int)'), partial(self.selectionChanged, axis=StatisticColumn.AXIS_Y))
        self.connect(self.checkbox_grid, SIGNAL('stateChanged(int)'), self.toggleGrid)
        self.connect(self.checkbox_legend, SIGNAL('stateChanged(int)'), self.toggleLegend)
        self.connect(self.save_plot, SIGNAL('clicked(bool)'), self.savePlot)
        self.connect(self.x_start, SIGNAL('textEdited(const QString&)'), partial(self.scaleChanged, start=True))
        self.connect(self.x_end, SIGNAL('textEdited(const QString&)'), partial(self.scaleChanged, start=False))
        
    def scaleChanged(self, string, start):
        if len(string) == 0:
            string = None
        if start is True:
            try:
                self.plot(x_start=string)
            except ValueError:
                # the intermediate user inserted string 
                # might not be a float after all
                pass
        else:
            try: 
                self.plot(x_end=string)
            except ValueError:
                pass
        
    def addColumn(self, statistic_column_id, description, axis):
        index = super(PlotController, self).addColumn(description, axis)
        if not index in self.__index2id.keys():
            self.__index2id[index] = {}
        self.__index2id[index][axis] = statistic_column_id
        
    def getLabel(self, id):
        unit = self.__columns[id]['unit']
        desc = self.__columns[id]['description']
        return "%s (%s)" % (desc, unit)
    
    def getData(self, simulation_id, statistic_column_id):
        data = self.__simulation_statistic_column.data(simulation_id, statistic_column_id)
        return data
    
    def __getSlice(self, data, start=None, end=None):
        """Return start and end index to catch all x in start<=x<end."""
        if start is None and end is None:
            return (0, -1)
        
        index_start = 0
        index_end = -1
        
        if start is not None:
            start = float(start)
        if end is not None:
            end = float(end)
            index_end = 0
            
        for item in data:
            item = float(item)
            if start is not None and item<start:
                index_start += 1
            if end is not None and item<end:
                index_end += 1
                
        return (index_start, index_end)
    
    def __calcNewScale(self, simulation_id, old_x_id, new_x_id, x_start, x_end):
        x_old = self.getData(simulation_id, old_x_id)
        x_new = self.getData(simulation_id, new_x_id)
        slice = self.__getSlice(x_old, x_start, x_end)
        return (x_new[slice[0]], x_new[slice[1]])
            
        
    def plot(self, x_id=None, y_id=None, grid=None, legend=None, x_start='unset', x_end='unset'):
        slice = None
        if x_id is None:
            x_id = self.__current_plot[StatisticColumn.AXIS_X]
        else:
            old_x_id = self.__current_plot[StatisticColumn.AXIS_X]
            if old_x_id is not None:
                scale = self.__calcNewScale(self.__simulations[0], old_x_id, x_id, self.__x_start, self.__x_end)
                self.x_start.setText(scale[0])
                self.x_end.setText(scale[1])
                x_start = scale[0]
                x_end = scale[1]
        if y_id is None:
            y_id = self.__current_plot[StatisticColumn.AXIS_Y]
        if grid is None:
            grid = self.__grid_state
        if legend is None:
            legend = self.__legend_state
        # unset instead of None intentionally!
        if x_start == 'unset':
            x_start = self.__x_start
        if x_end == 'unset':
            x_end = self.__x_end
            
        simulation = Simulation() 
        self.canvas.clear()
        self.canvas.setLabel(x=self.getLabel(x_id), y=self.getLabel(y_id))
        for simulation_id in self.__simulations:
            x = self.getData(simulation_id, x_id)
            y = self.getData(simulation_id, y_id)
            if x is None or y is None:
                continue
            if slice is None:
                slice = self.__getSlice(x, x_start, x_end)           
            x = x[slice[0]:slice[1]]
            y = y[slice[0]:slice[1]]              
            label = simulation_id
            note = simulation.note(simulation_id)
            if len(note['note_title']) > 0:
                label += ' - %s' % note['note_title']
            self.canvas.addPlot(x, y, label=label)
        self.canvas.grid(grid)
        if legend:
            self.canvas.legend()
        self.canvas.draw()
        pyplot.show()
        self.__current_plot = {StatisticColumn.AXIS_X: x_id, StatisticColumn.AXIS_Y: y_id}
        self.__grid_state = grid
        self.__legend_state = legend
        self.__x_end = x_end
        self.__x_start = x_start
        
    def selectionChanged(self, combo_box_index, axis):
        id = self.__index2id[combo_box_index][axis]
        if axis == StatisticColumn.AXIS_X: 
            self.plot(x_id=id)
        else:
            self.plot(y_id=id)
            
    def toggleGrid(self, state):
        """
        TODO: must be possible to toggle the grid without a complete re-plot.
        """
        if state == 2:
            self.plot(grid=True)
        else:
            self.plot(grid=False)
            
    def toggleLegend(self, state):
        """
        TODO: must be possible to toggle the legend without a complete re-plot
        """
        if state == 2:
            self.plot(legend=True)
        else:
            self.plot(legend=False)
    
    def savePlot(self, state):
        simulation = Simulation() 
        x_id = self.__current_plot['x']
        y_id = self.__current_plot['y']
        x_name = self.getLabel(x_id)
        y_name = self.getLabel(y_id)
        for simulation_id in self.__simulations:
            x = self.getData(simulation_id, x_id)
            y = self.getData(simulation_id, y_id)
            if x is None or y is None:
                continue
            label = simulation_id
            note = simulation.note(simulation_id)
            if len(note['note_title']) > 0:
                label += ' - %s' % note['note_title']
        label += " - " + x_name + " vs. " + y_name
        pyplot.savefig(label + ".png")
            
class VisibleColumnPickerController(QDialog, VisibleColumnPickerView):
    
    def __init__(self, accelerator_id, parent=None):
        QDialog.__init__(self, parent)
        VisibleColumnPickerView.__init__(self)
        
        self.__parent = parent
        self.__accelerator_id = accelerator_id
        self.__group2parameter = {}     # {group_id: [parameter_id, ...], ...}
        self.__group_visible = {}       # {group_id: state, ...}
        self.__column_visible = {}     # {parameter_id: state, ...} not influenced by the group check state!
        
        self.setupUi()
        self.readSettings()
        if parent is not None:
            self.resize(parent.width()*0.8, parent.height()*0.8)
            
        accelerator = Accelerator()
        parameters = accelerator.parameters(accelerator_id)
        
        self.setWindowTitle('OPAL Pilot > %s > Column Picker' % accelerator.alias(accelerator_id))
        
        group = self.addGroup(0, 'General', self.isGroupVisible(0))
        self.connect(group, SIGNAL('toggled(bool)'), partial(self.toggleGroup, group_id=0))
        general_columns = {
            -7: 'Date',
            -6: 'Username',
            -5: 'Note',
            -4: 'State',
            -3: 'ID',
            -2: 'Process',
            -1: 'Server',
             0: 'Cores'
        }
        self.__group2parameter[0] = []
        for parameter_id in sorted(general_columns.iterkeys()):
            label = general_columns[parameter_id]
            column = self.addColumn(0, label, self.isColumnVisible(parameter_id))
            self.connect(column, SIGNAL('stateChanged(int)'), partial(self.toggleColumn, parameter_id=parameter_id))
            self.__group2parameter[0].append(parameter_id)
        
        if parameters is not None:
            groups = []
            for parameter in parameters:
                gid = int(parameter['group_id'])
                pid = int(parameter['id'])
                
                if not gid in self.__group2parameter.keys():
                    group = self.addGroup(gid, parameter['group_alias'], self.isGroupVisible(gid))
                    self.__group2parameter[gid] = []
                    self.connect(group, SIGNAL('toggled(bool)'), partial(self.toggleGroup, group_id=gid))
                    
                column = self.addColumn(gid, parameter['label'], self.isColumnVisible(pid))
                self.connect(column, SIGNAL('stateChanged(int)'), partial(self.toggleColumn, parameter_id=pid))
                self.__group2parameter[gid].append(pid)
                
        self.connect(self.buttons, SIGNAL('accepted()'), self.accept)              
                
        self.show()
        
    def toggleGroup(self, checked, group_id):
        self.__group_visible[group_id] = checked
        if checked is False:
            for parameter_id in self.__group2parameter[group_id]:
                self.__parent.toggleColumn(parameter_id, column_visible=False)
        else:
            for parameter_id in self.__group2parameter[group_id]:
                if parameter_id not in self.__column_visible.keys():
                    state = True
                else:
                    state = self.__column_visible[parameter_id]
                self.__parent.toggleColumn(parameter_id, state)
                         
    def toggleColumn(self, state, parameter_id):
        if parameter_id not in self.__column_visible.keys():
            self.__column_visible[parameter_id] = False
        else:
            self.__column_visible[parameter_id] = not self.__column_visible[parameter_id]
        self.__parent.toggleColumn(parameter_id)
        
    def isGroupVisible(self, group_id):
        if group_id not in self.__group_visible.keys():
            return True
        return self.__group_visible[group_id]
    
    def isColumnVisible(self, parameter_id):
        if parameter_id not in self.__column_visible.keys():
            return True
        return self.__column_visible[parameter_id]
    
    def accept(self):
        self.close()
        
    def closeEvent(self, event):
        self.saveSettings()
        event.accept()
        
    def saveSettings(self):
        settings = QSettings()
        settings.beginGroup('AcceleratorVisibleColumns')
        
        settings.beginWriteArray(self.__accelerator_id + '_groups')
        i = 0
        for group_id in self.__group_visible.keys():
            settings.setArrayIndex(i)
            settings.setValue('group_id', group_id)
            settings.setValue('state', self.__group_visible[group_id])
            i += 1
        settings.endArray()
        
        settings.beginWriteArray(self.__accelerator_id + '_columns')
        i = 0
        for parameter_id in self.__column_visible.keys():
            settings.setArrayIndex(i)
            settings.setValue('parameter_id', parameter_id)
            settings.setValue('state', self.__column_visible[parameter_id])
            i += 1
        settings.endArray()

        settings.endGroup()
        
    def readSettings(self):
        settings = QSettings()
        settings.beginGroup('AcceleratorVisibleColumns')
        size = settings.beginReadArray(self.__accelerator_id + '_groups')
        for i in range(size):
            settings.setArrayIndex(i)
            group_id = settings.value('group_id').toString()
            state = settings.value('state').toBool()
            self.__group_visible[int(group_id)] = state
        settings.endArray()
        size = settings.beginReadArray(self.__accelerator_id + '_columns')
        for i in range(size):
            settings.setArrayIndex(i)
            parameter_id = settings.value('parameter_id').toString()
            state = settings.value('state').toBool()
            self.__column_visible[int(parameter_id)] = state        
        settings.endArray()        
        settings.endGroup()
        

class InputMaskController(QDialog, InputMaskView):
    
    """
    Render generic input form.
    
    The form is based on the database tables:
        accelerator_parameters
        parameters
        server
        ui_group
        ui_unit
        
    """
    
    """
    TODO: load "default" values based on simulation_id
    TODO: make dialog blocking
    TODO: validate values from the input fields 
    """        
    def __init__(self, accelerator_id, parent=None, simulation_id=None):      
        QDialog.__init__(self, parent)
        InputMaskView.__init__(self)
        
        self.__parent = parent
        self.__parameter_map = {}   # Mapping parameter_id to QLineEdit
        self.__server_map = {}      # Mapping QComboBox index to server_id
        
        self.setupUi()
        
        if parent is not None:
            self.resize(parent.width()*0.8, parent.height()*0.8)         
        
        simulation = Simulation()
        accelerator = Accelerator()
        server = Server()
        parameters = accelerator.parameters(accelerator_id, simulation_id)
        servers = server.names()
        self.__accelerator_id = accelerator_id
        
        if simulation_id is None:
            base = 'Default'
        else:
            base = simulation_id
        self.setWindowTitle('OPAL Pilot > %s > New Simulation (Based on: %s)' % (accelerator.alias(accelerator_id), base))
                
        
        groups = []
        for parameter in parameters:
            gid = parameter['group_id']
            pid = parameter['id']
            
            if not gid in groups:
                groups.append(gid)
                self.addGroup(gid, parameter['group_alias'])
                
            line_edit = self.addParameter(gid, pid, parameter['default_value'], parameter['label'], parameter['unit'], parameter['tooltip'])
            self.__parameter_map[pid] = line_edit
            
        for server in servers:
            index = self.addServer(server['id'], server['ui_alias'])
            self.__server_map[index] = server['id']
            
        if simulation_id is not None:
            """ TODO: pick correct server """
            note = simulation.note(simulation_id)
            server = simulation.server(simulation_id)
            self.setDefault(cores=server['server_cores'], note_title=note['note_title'], note=note['note'])
                                            
        self.connect(self.buttons, SIGNAL('accepted()'), self.accept)
        self.connect(self.buttons, SIGNAL('rejected()'), self.reject)
        self.show()
          
    def accept(self):   
        """ Launch simulation and return simulation_id """
        parameters = {}
        for pid in self.__parameter_map.iterkeys():
            """ TODO: check if empty or even validate properly """
            parameters[pid] = self.__parameter_map[pid].text()
            
        server_index = self.server()
        server_id = self.__server_map[server_index]
        server_cores = self.serverCores()
        
        note_title, note = self.note()
                    
        simulation = Simulation()
        simulation_id = simulation.launch(self.__accelerator_id, server_id, server_cores, note_title, note, parameters)
        if self.__parent is not None:
            self.__parent.newRow(simulation_id)        
        
        self.close()
