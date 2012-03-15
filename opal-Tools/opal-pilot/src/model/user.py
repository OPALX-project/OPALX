import os
import pwd

from . import Model

class User(Model):
    
    def __init__(self):
        Model.__init__(self)
        
    def id(self):
        username = self.__getUsername()
        
        query = """ SELECT id FROM %s WHERE username = '%s' """ % (self.table, username)
        id = self.conn.field(query)
        
        if id is None:
            id = self.conn.insert('user', ['username'], [(username)])
        
        return id
        
    def __getUsername(self):
        return pwd.getpwuid(os.getuid()).pw_name