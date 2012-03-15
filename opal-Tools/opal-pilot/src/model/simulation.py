import re
import string
import tempfile
import warnings

from . import Model, CachingModel
from .accelerator import Accelerator
from .server import Server
from .user import User

class SimulationStatisticColumn(CachingModel):
    
    VALUE_SEPERATOR = ';'
    
    def __init__(self):
        CachingModel.__init__(self)
        self.__insert_stack = {}    # {simulation_statistic_column_id: [value, value, ...], ...}
    
    def id(self, simulation_id, statistic_column_id, use_cache=True, create_if_not_found=True):
        """Retrieve and return the simulation statistic column id.
        
        None will be returned if no id can be found and create_if_not_found
        is not True.
        
        """
        if use_cache is True:
            id = self.cache.find("%s,%s" % (simulation_id, statistic_column_id))
            if id is not None:
                return id
            
        query = """
        SELECT id 
        FROM %s 
        WHERE simulation_id = '%s' 
        AND statistic_column_id = '%s'
        """ % (self.table, simulation_id, statistic_column_id)
        id = self.conn.field(query)
        
        if id is None and create_if_not_found is True:
            id = self.conn.insert(self.table, ['simulation_id', 'statistic_column_id'], [(simulation_id, statistic_column_id)])
            
        if id is not None:
            self.cache.put("%s,%s" % (simulation_id, statistic_column_id), id)
            
        return id
    
    def data(self, simulation_id, statistic_column_id, use_cache=True):
        simulation_statistic_column_id = self.id(simulation_id, statistic_column_id, create_if_not_found=False)
        if simulation_statistic_column_id is None:
            return None
        
        if use_cache is True:
            data = self.cache.find('data:%s' % simulation_statistic_column_id)
            if data is not None:
                return data
        
        query = """SELECT `values` FROM %s WHERE id = %s""" % (self.table, simulation_statistic_column_id)
        data = self.conn.field(query)
        data = data.split(self.VALUE_SEPERATOR)
        self.cache.put('data:%s' % simulation_statistic_column_id, data)
        return data
    
    def columnIds(self, simulation_ids):
        query = """
        SELECT DISTINCT statistic_column_id
        FROM %s
        WHERE simulation_id IN (%s)
        """ % (self.table, ','.join(simulation_ids))
        return self.conn.rs(query, how=self.conn.RETURN_TUPLE)
    
    def valuePrepare(self, simulation_statistic_column_id, value):
        """Prepare a single value (one item from a column) for insertion."""
        if not simulation_statistic_column_id in self.__insert_stack.keys():
            self.__insert_stack[simulation_statistic_column_id] = []
        self.__insert_stack[simulation_statistic_column_id].append(value)
        
    def valueInsert(self, simulation_statistic_column_id):
        """Insert all values from the insertion stack to the database.
        
        Values are stacked up with the valuePrepare method.
        
        """
        if not simulation_statistic_column_id in self.__insert_stack.keys():
            warnings.warn('Insert was issued without data to insert.', Warning, stacklevel=2)
            return
        
        stack = self.__insert_stack[simulation_statistic_column_id]
        self.conn.execute("UPDATE %s SET `values` = NULL WHERE id = %s" % (self.table, simulation_statistic_column_id))
        self.conn.blobAppend(self.table, 'values', stack, simulation_statistic_column_id, self.VALUE_SEPERATOR)
        del self.__insert_stack[simulation_statistic_column_id]
        

class Simulation(Model):
    
    STATE_QUEUED = '1 - Queued'
    STATE_RUNNING = '2 - Running'
    STATE_POSTPROCESSING_QUEUE = '3 - Queued for Postprocessing'
    STATE_POSTPROCESSING = '4 - Postprocessing'
    STATE_FINISHED = '5 - Finished'
    STATE_ABORTED = '6 - Aborted'
    
    States = (STATE_QUEUED, STATE_RUNNING, STATE_POSTPROCESSING_QUEUE, STATE_POSTPROCESSING, STATE_FINISHED, STATE_ABORTED)
    StatesUnfinished = (STATE_QUEUED, STATE_RUNNING)
    
    __TEMPLATE_VARIABLE_IDENTIFIER = '_%s_'
    __TEMP_DIR = 'tmp/'
    
    __PROCESS_ID_PATTERN = 'Your job ([0-9]*) \(".*"\) has been submitted'
    
    def __init__(self):
        Model.__init__(self)
    
    def setState(self, process_id, state):
        """ Update the state of a simulation. """
        if state not in self.States:
            raise Exception('Invalid state given')
        query = """ UPDATE %s SET state = '%s' WHERE process_id = %s """ % (self.table, state, process_id)
        self.conn.execute(query)
        
    def getState(self, simulation_id):
        query = """ SELECT state FROM %s WHERE id = %s """ % (self.table, simulation_id)
        return self.conn.field(query)
    
    def unfinished(self):
        """ 
        Determine unfinished simulations and return dictionary.
        
        Unfinished in the means of: Not completely simulated by OPAL yet.  
        
        The dictionary will be in the form:
        {process_id: {'id': id, 'accelerator_id': accelerator_id, 'server_id': server_id}, ...} 
        
        """
        query = """
        SELECT id, server_id, accelerator_id, process_id
        FROM %s
        WHERE deleted = 0
        AND state IN ('%s')
        """ % (self.table, "','".join(self.StatesUnfinished))
        result = self.conn.rs(query, how=self.conn.RETURN_TUPLE)
        
        if result is None:
            return {}
        
        output = {}
        for item in result:
            simulation_id, server_id, accelerator_id, process_id = item
            output[process_id] = {
                'server_id': server_id,
                'accelerator_id': accelerator_id,
                'id': simulation_id
            }
            
        return output
    
    def nextToParse(self):
        """
        Determine the next simulation that needs stat file parsing.
        """
        query = """ SELECT id, server_id, process_id FROM %s WHERE state = '%s' ORDER BY id LIMIT 0,1 """ % (self.table, self.STATE_POSTPROCESSING_QUEUE)
        return self.conn.row(query, how=self.conn.RETURN_DICT) 
        
                      
    def ids(self, accelerator_id, last_updated_ms=None):
        """ Determine simulation_ids for the given accelerator and return a list. """
        if last_updated_ms is None:
            query = """SELECT id FROM %s WHERE deleted = 0 AND accelerator_id = %s """ % (self.table, accelerator_id)
        else:
            query = """SELECT id FROM %s WHERE deleted = 0 AND accelerator_id = %s
                    AND updated > NOW() - INTERVAL %s SECOND""" % (self.table, accelerator_id, last_updated_ms/1000)
        simulation_ids = []
        result = self.conn.rs(query, how=self.conn.RETURN_TUPLE)
        if result is None:
            return simulation_ids
        for item in result:
            id = item[0]
            simulation_ids.append(id)
        return simulation_ids
    
    def note(self, simulation_id):
        """Return the note_title and note of a simulation"""
        query = """SELECT note, note_title FROM %s WHERE id = %s""" % (self.table, simulation_id)
        return self.conn.row(query, how=self.conn.RETURN_DICT)
        
    def server(self, simulation_id):
        """Return the server related information of a simulation"""
        query = """SELECT server_id, server_cores FROM %s WHERE id = %s""" % (self.table, simulation_id)
        return self.conn.row(query, how=self.conn.RETURN_DICT)
        
    
    def values(self, simulation_ids):
        """ Determine the input values for one or more simulation_ids. """ 
        output = {} # mapping simulation_id -> values
        if simulation_ids.__class__ is long:
            simulation_ids = [simulation_ids]
            warnings.warn('Expecting simulation_ids to be a list, but got long', category=SyntaxWarning, stacklevel=2)
        simulation_ids = [str(x) for x in simulation_ids]
        
        query = """
        SELECT %(table)s.id AS id, server.ui_alias AS server, server.id AS server_id, process_id, DATE_FORMAT(time, '%%d.%%m.%%Y %%H:%%i') AS date, state, server_cores, user.username AS username, note_title, note
        FROM %(table)s, server, user
        WHERE %(table)s.id IN (%(ids)s)
        AND %(table)s.server_id = server.id
        AND %(table)s.user_id = user.id
        ORDER BY id DESC
        """ % {'table': self.table, 'ids': ','.join(simulation_ids)}
        simulations = self.conn.rs(query)
        
        if simulations is None:
            return output
          
        for simulation in simulations:
            output[simulation['id']] = {
                'server': simulation['server'],
                'server_id': simulation['server_id'],
                'process_id': simulation['process_id'],
                'date': simulation['date'],
                'state': simulation['state'],
                'server_cores': simulation['server_cores'],
                'username': simulation['username'],
                'note': simulation['note'],
                'note_title': simulation['note_title'],
                'values': {}
            } 
        
        query = """
        SELECT simulation_id, parameter_id, value
        FROM simulation_parameter
        WHERE simulation_id IN (%s)
        """ % ','.join(simulation_ids)
        parameters = self.conn.rs(query, how=self.conn.RETURN_TUPLE)
        
        for parameter in parameters:
            simulation_id, parameter_id, value = parameter
            output[simulation_id]['values'][parameter_id] = value       

        return output
        
    def launch(self, accelerator_id, server_id, server_cores, note_title, note, parameter_values):
        """
        Launch new simulation and return simulation_id
        
        Keyword arguments:
        accelerator_id
        server_id
        server_cores -- CPU cores to use on the Server
        parameter_values -- Dictionary in the form {id: values, ...}
        
        """
        self.__validateId('accelerator', accelerator_id, 'Invalid accelerator (ID: %s)')
        self.__validateId('server', server_id, 'Invalid server (ID: %s)')
        self.__validateId('parameter', parameter_values.keys(), 'Invalid parameters')
        
        accelerator = Accelerator()
        server = Server()
        user = User()
        ssh = server.connection(server_id)

        user_id = user.id()
        print "--->", server_id, accelerator_id, server_cores, parameter_values, user_id        
        simulation_id = self.__insertToDB(server_id, accelerator_id, server_cores, note_title, note, parameter_values, user_id)
        
        data_dir = '%s/' % server.path(server_id)
        simulation_dir = data_dir + '%s/' % simulation_id
        filenames = accelerator.filenames(accelerator_id)
        
        ssh.execute('mkdir %s' % simulation_dir)
        # Write infile       
        content = self.__parseInfile(accelerator_id, parameter_values)
        self.__put(ssh, simulation_dir + 'opalgui.in', content)
        # Write phasefile
        content = accelerator.fileContent(accelerator_id, 'opalgui.phases')
        self.__put(ssh, simulation_dir + 'opalgui.phases', content)
        # Upload and link files
        for filename in filenames:
            # Upload only files that are not present on the server
            if not ssh.isfile(data_dir + filename):
                content = accelerator.fileContent(accelerator_id, filename)
                self.__put(ssh, data_dir + filename, content)
            ssh.link('../' + filename, simulation_dir)

        # Start simulation      
        output = ssh.execute('cd %(simulation_dir)s && qsub -q all.q -pe mpi %(cores)s %(filename)s' % {'simulation_dir': simulation_dir, 'cores': server_cores, 'filename': 'opalgui.sge'})
        print output[-1]
        match = re.match(self.__PROCESS_ID_PATTERN, output[-1])
        if match is not None:
            process_id = match.group(1)
        self.conn.execute('UPDATE %(table)s SET process_id = %(process_id)s WHERE id = %(simulation_id)s' % {'table': self.table, 'process_id': process_id, 'simulation_id': simulation_id})
        
        return simulation_id
    
    def __parseInfile(self, accelerator_id, values):
        accelerator = Accelerator()
        rs = accelerator.parameters(accelerator_id)
        parameters = {}
        for item in rs:
            parameters[item['id']] = item['parameter']

        infile = accelerator.fileContent(accelerator_id, 'opalgui.in')                            
        for id in values.iterkeys():
            infile = infile.replace(self.__TEMPLATE_VARIABLE_IDENTIFIER % parameters[id], str(values[id]))
            
        return infile
           
    def __validateId(self, database_table, id, error_code):
        if not self.conn.hasId(database_table, id):
            print 'todo debug:', id.__class__
            if id.__class__ is str or id.__class__ is int:
                raise SimulationError(error_code % str(id))
            else:
                raise SimulationError(error_code)
            
    def __insertToDB(self, server_id, accelerator_id, server_cores, note_title, note, parameter_values, user_id):
        columns = ['server_id', 'accelerator_id', 'time', 'state', 'server_cores', 'user_id', 'note_title', 'note']
        values = [(server_id, accelerator_id, self.conn.now(), self.STATE_QUEUED, server_cores, user_id, note_title, note)]
        simulation_id = self.conn.insert(self.table, columns, values)
        
        columns = ['simulation_id', 'parameter_id', 'value']
        values = []
        for id in parameter_values.iterkeys():
            values.append((simulation_id, id, parameter_values[id]))
        self.conn.insert('simulation_parameter', columns, values)
        return simulation_id
    
    def __put(self, ssh, remote_path, content):
        temp = tempfile.NamedTemporaryFile('wb')
        temp.write(content)
        temp.flush()
        ssh.put(temp.name, remote_path)
        temp.close()    
                    
class SimulationError():
    
    def __init__(self, msg):
        self.msg = msg