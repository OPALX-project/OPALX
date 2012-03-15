from . import Model
from util.decorator import deprecated

class Accelerator(Model):
    
    def __init__(self):
        Model.__init__(self)
            
    def names(self):
        query = """ SELECT id, ui_alias FROM %s ORDER BY id """ % (self.table)
        return self.conn.rs(query)
    
    def alias(self, accelerator_id):
        query = """ SELECT ui_alias FROM %s WHERE id = %s """ % (self.table, accelerator_id)
        return self.conn.field(query)
        
    def parameters(self, accelerator_id, simulation_id=None):
        self.__checkId(accelerator_id)        
        
        if simulation_id is None:        
            query = """
            SELECT id, parameter, label, group_id, group_alias, unit, tooltip, default_value 
            FROM view_accelerator_parameter
            WHERE accelerator_id = %(acc_id)s
            """ % {'acc_id': accelerator_id}
        else:
            query = """
            SELECT id, parameter, label, group_id, group_alias, unit, tooltip, (
                SELECT value 
                FROM simulation_parameter
                WHERE simulation_id = %(sim_id)s
                AND simulation_parameter.parameter_id = view_accelerator_parameter.id 
            ) AS default_value
            FROM view_accelerator_parameter
            WHERE accelerator_id = %(acc_id)s
            """ % {'acc_id': accelerator_id, 'sim_id': simulation_id}
        return self.conn.rs(query, how=self.conn.RETURN_DICT)
    
    def __file(self, accelerator_id, type):
        self.__checkId(accelerator_id)
        
        query = """
        SELECT %s
        FROM %s
        WHERE id = %s
        """ % (type, self.table, accelerator_id)
        return self.conn.field(query)
    
    def filenames(self, accelerator_id):
        self.__checkId(accelerator_id)
        
        filenames = []
        
        query = """
        SELECT file.filename as filename
        FROM file, accelerator_file
        WHERE accelerator_file.accelerator_id = %s 
        AND accelerator_file.file_id = file.id        
        """ % accelerator_id
        rs = self.conn.rs(query, how=self.conn.RETURN_TUPLE)
        
        if rs is not None:        
            for item in rs:
                filenames.append(item[0])
                
        return filenames
    
    def fileContent(self, accelerator_id, filename):
        if filename == 'opalgui.in':
            return self.__file(accelerator_id, 'input_file')
        elif filename == 'opalgui.phases':
            return self.__file(accelerator_id, 'phase_file')

        query = """
        SELECT file.content as content
        FROM file, accelerator_file
        WHERE accelerator_file.accelerator_id = %s
        AND file.filename = '%s'
        AND accelerator_file.file_id = file.id
        """ % (accelerator_id, filename)
        
        content = self.conn.field(query)
        return content
            
        
    @deprecated
    def files(self, accelerator_id):
        self.__checkId(accelerator_id)
        
        query = """
        SELECT input_file, phase_file
        FROM accelerator
        WHERE id = %s
        """ % accelerator_id
        input_file, phase_file = self.conn.row(query, self.conn.RETURN_TUPLE)
        
        files = {}
        files['opalgui.in'] = input_file
        files['opalgui.phases'] = phase_file
        
        query = """ 
        SELECT file.content as content, file.filename as filename
        FROM file, accelerator_file
        WHERE accelerator_file.accelerator_id = %s AND
        accelerator_file.file_id = file.id
        """ % accelerator_id
        rs = self.conn.rs(query)
        
        for item in rs:
            files[item['filename']] = item['content']   
        return files        
        
    def __checkId(self, accelerator_id):
        if not self.conn.hasId(self.table, accelerator_id):
            raise Exception('Unknown Accelerator (id: %s)' % accelerator_id)
        return True
            
     
        
        
        