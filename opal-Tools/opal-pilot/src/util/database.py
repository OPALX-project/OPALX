import _mysql

from util.decorator import hotshotit

class DatabaseConnection():
    
    """
    A layer for the _mysql interface to make working with it more convenient.
    
    The connection is made persistent by making the class a singleton.
    Interaction with the database is enable through the following methods:
    -- Inserting
    insert - Insert one or more rows into the database table.
    """
    """TODO: Read connection parameters from config file."""
    
    """TODO: Persistence by singleton does not work together with threads.""" 
    # __metaclass__ = Singleton
    
    ROWS_ALL = 0
    RETURN_TUPLE = 0
    RETURN_DICT = 1
    
    __DEBUG_QUERIES = False
    
    def __init__(self):
        self.db = {
            'user': 'opaldb_rw',
            'pass': 'is4opaldb',
            'db': 'OPALDB', 
            'port': 3306,
            'host': 'pythia-alpha.psi.ch'
        }
        self.conn = None 
    
    def connect(self):
        if self.conn is None:
            self.conn = _mysql.connect(host=self.db['host'], port=self.db['port'], db=self.db['db'], user=self.db['user'], passwd=self.db['pass'])             
        
    def blobAppend(self, table, field, values, id, seperator=None):
        query = """ UPDATE `%(table)s` SET `%(field)s` = %(concat_function)s WHERE id = %(id)s """
        mappings = {'table': table, 'field': field, 'id': id} 
        if seperator is None:
            mappings['concat_function'] = """CONCAT(`%s`, '%s')""" % (field, ''.join(values))
        else:
            mappings['concat_function'] = """CONCAT_WS('%s', `%s`, '%s')""" % (seperator, field, seperator.join(values))
        self.execute(query % mappings) 
        
    def execute(self, query):
        self.connect()
        self.conn.query(query)
      
    def insert(self, table, fields, values, max_rows=None):
        """
        Insert new row(s) into table and return insert id.
        
        Inserting is very slow. Combine as much rows per insert as possible.
        
        Keyword arguments:
        table -- the database table
        fields -- list of database field names
        values -- list of tuples containing the values. i.e.: [('first', 2, 4), ...]
        max_rows -- number of rows to insert at a time
         
        """
        self.connect()
        print table
        print fields
        print values
        print "--"
        
        # Split to max_rows chunks
        if isinstance(max_rows, int) and max_rows > 0:
            while len(values) > 0:
                self.insert(table, fields, values[:max_rows])
                values = values[max_rows:]

        size = len(fields)
        if not isinstance(values, list):
            raise DatabaseConnectionError('Values parameter is expected to be a list of tuples. Got %s instead' % values.__class__)
        if not isinstance(values[0], tuple) and size > 1:
            raise DatabaseConnectionError('Expecting tuples in values list, but got %s' % values[0].__class__)
        if size == 1:
            values = ['(' + self.__prepareValue(value) + ')' for value in values]
        else:
            for i in range(len(values)):
                if len(values[i]) != size:
                    raise DatabaseConnectionError('Tuple in values list does not match fields length: %s' % values[i])
                values[i] = [self.__prepareValue(value) for value in values[i]]
                values[i] = "(" + ",".join(values[i]) + ")"

        mappings = {'table': table, 'fields': ','.join(fields), 'values': ','.join(values)}
        query = 'INSERT INTO %(table)s (%(fields)s) VALUES %(values)s' % mappings
        self.conn.query(query) 
        return self.conn.insert_id()  
            
    def rs(self, query, maxrows=ROWS_ALL, how=RETURN_DICT):
        """ 
        Submit SELECT statement and return results.
        
        Keyword arguments:
        query -- SELECT statement
        maxrows -- rows to return. default: ROWS_ALL
        how -- format of the return. default: RETURN_DICT.  
               RETURN_DICT returns dictionarys in a list [{'rowname': value, ...}, ...]
               RETURN_TUPLE returns tuples in a list [(value, ...), ...]
        
        """
        self.connect()
        self.conn.query(query)        
        result = self.conn.store_result().fetch_row(maxrows=maxrows, how=how)
        
        if self.__DEBUG_QUERIES:
            self.__explain(query)
                    
        if len(result) == 0:
            return None
        return result
    
    def field(self, query):
        """ Return the first field of the SELECT statement. """
        result = self.row(query, how=self.RETURN_TUPLE)
        if result is None or len(result) == 0:
            return None
        return result[0]

    def row(self, query, how=RETURN_DICT):
        """ Return the first row of the SELECT statement. """
        result = self.rs(query, maxrows=1, how=how)
        if result is None or len(result) == 0:
            return None
        return result[0]
    
    def now(self):
        """ Return result of the MySQL NOW() function. """
        query = 'SELECT NOW()'
        return self.field(query)                
    
    def hasId(self, table, id):
        """
        Return if one or more ids exist in a table.
        
        Keyword arguments:
        table -- database table
        id -- int or list of id(s)
        
        """
        self.connect()
        if id.__class__ is str:
            id = [id]
        elif id.__class__ is list:
            pass
        else:
            raise DatabaseConnectionError('Wrong format for parameter id: %s' % id.__class__)
            
        ids = ','.join(id)
        self.conn.query("SELECT id FROM %(table)s WHERE id IN (%(id)s)" % {'table': table, 'id': ids})
        if self.conn.store_result().num_rows() < len(id):
            return False
        return True
    
    def escape_string(self, string):
        self.connect()
        return self.conn.escape_string(string)
    
    def __prepareValue(self, value):
        if value.__class__ is int:
            return str(value)
        else:
            return '"%s"' % value
        
    def __explain(self, query):
        self.__DEBUG_QUERIES = False
        explain = self.row("EXPLAIN " + query)
        self.__DEBUG_QUERIES = True
        print  
        print "-- EXPLAIN QUERY ----"
        print query
        for key in explain.iterkeys():
            print "   %s: %s" % (key, explain[key])         

    
class DatabaseConnectionError():
    
    def __init__(self, msg):
        self.msg = msg
        
    def __str__(self):
        return repr(self.msg)