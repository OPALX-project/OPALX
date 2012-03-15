import re

from . import Model
from . import Cache
from .simulation import SimulationStatisticColumn 

class StatisticType(Model):
    
    def __init__(self):
        Model.__init__(self)
        self.__cache = Cache()
        
    def id(self, type, use_cache=True, create_if_not_found=True):
        """Retrieve and return the id.
        
        None will be returned if no id can be found and create_if_not_found
        is not True.
        
        """
        if use_cache is True:
            id = self.__cache.find(type)
            if id is not None:
                return id
            
        query = """SELECT id FROM %s WHERE type = '%s'""" % (self.table, type)
        print query
        id = self.conn.field(query)
        
        if id is None and create_if_not_found is True:
            id = self.conn.insert(self.table, ['type'], [(type)])
            
        if id is not None:
            self.__cache.put(type, id)
            
        return id
    
class StatisticUnit(Model):
    
    def __init__(self):
        Model.__init__(self)
        self.__cache = Cache()
        
    def id(self, unit, use_cache=True, create_if_not_found=True):
        """Retrieve and return the unit id.
        
        None will be returned if no id can be found and create_if_not_found
        is not True.
        
        """
        if use_cache is True:
            id = self.__cache.find(unit)
            if id is not None:
                return id
            
        query = """SELECT id FROM %s WHERE unit = '%s'""" % (self.table, unit)
        id = self.conn.field(query)
        
        if id is None and create_if_not_found is True:
            id = self.conn.insert('statistic_unit', ['unit'], [(unit)])
            
        if id is not None:
            self.__cache.put(unit, id)
            
        return id
    
class StatisticColumn(Model):
    
    AXIS_X = 'x'
    AXIS_Y = 'y'
    
    def __init__(self):
        Model.__init__(self)
        self.__cache = Cache()
        self.type = StatisticType()
        self.unit = StatisticUnit() 
        
    def id(self, name, type, unit, description, use_cache=True, create_if_not_found=True):
        """Retrieve and return the column id.
        
        None will be returned if no id can be found and create_if_not_found
        is not True.
        
        """
        type_id = self.type.id(type)
        unit_id = self.unit.id(unit)
        if use_cache is True:
            id = self.__cache.find("%s,%s,%s,%s" % (name, type_id, unit_id, description))
            if id is not None:
                return id
            
        query = """
        SELECT id 
        FROM %s
        WHERE name = '%s'
        AND statistic_type_id = %s
        AND statistic_unit_id = %s
        AND description = '%s'
        """ % (self.table, name, type_id, unit_id, description)
        id = self.conn.field(query)
        
        if id is None and create_if_not_found is True:
            id = self.conn.insert('statistic_column', ['name', 'statistic_type_id', 'statistic_unit_id', 'description'], [(name, type_id, unit_id, description)])        
        
        if id is not None:
            self.__cache.put("%s,%s,%s,%s" % (name, type_id, unit_id, description), id)
            
        return id
    
    def columns(self, statistic_column_ids):
        query = """
        SELECT %(table)s.id AS id, name, type, unit, description, ui_axis
        FROM %(table)s, statistic_type, statistic_unit
        WHERE %(table)s.id IN (%(ids)s)
        AND %(table)s.statistic_type_id = statistic_type.id
        AND %(table)s.statistic_unit_id = statistic_unit.id
        ORDER BY description
        """ % {'table': self.table, 'ids': ','.join(statistic_column_ids)}
        return self.conn.rs(query, how=self.conn.RETURN_DICT)

              
class Statistic(Model):
    
    def __init__(self):
        Model.__init__(self)
        self.statistic_column = StatisticColumn()
        self.simulation_statistic_column = SimulationStatisticColumn()
        
    def parseStatFile(self, simulation_id, content):
        regexp = {
            'settings': '(&.*&end)',
            'columns': '&column name=(.*?)[,\s]+type=(.*?)[,\s]+units=(.*?)[,\s]+description="(.*?)"',
            'parameters': '&parameter name=(.*?)[,\s]+',
            'data': '.*&end(.*)'
        }
        
        statfile = {}
        statfile['settings'] = re.search(regexp['settings'], content, re.DOTALL).group(1)
        statfile['columns'] = re.findall(regexp['columns'], statfile['settings'], re.DOTALL)
        statfile['parameters'] = re.findall(regexp['parameters'], statfile['settings'], re.DOTALL)
        statfile['content'] = re.search(regexp['data'], content, re.DOTALL).group(1)
        statfile['content'] = statfile['content'].split('\n')
        start = len(statfile['parameters']) + 1
        end = len(statfile['content']) - 1
        statfile['content'] = statfile['content'][start:end]
        
        row = statfile['content'][0]
        items = row.split('\t')
        num_cols = len(items) - 1
        num_rows = len(statfile['content'])
        matrix = [[0 for i in range(num_cols)] for j in range(num_rows)]
    
        for i in range(num_rows):
            items = [item.strip() for item in statfile['content'][i].split('\t')]
            for j in range(num_cols):
                matrix[i][j] = items[j]
                
        # Select a single column from the matrix. 
        for j in range(num_cols):
            name, type, unit, desc = statfile['columns'][j]
            statistic_column_id = self.statistic_column.id(name, type, unit, desc)
            simulation_statistic_column_id = self.simulation_statistic_column.id(simulation_id, statistic_column_id)
            # Select a single row from the column
            for i in range(num_rows):
                self.simulation_statistic_column.valuePrepare(simulation_statistic_column_id, matrix[i][j])
            self.simulation_statistic_column.valueInsert(simulation_statistic_column_id)