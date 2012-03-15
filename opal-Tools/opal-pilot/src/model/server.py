from . import Model
from util.ssh import SSHConnection

class Server(Model):
    
    def names(self):
        """ Retrieve servers and return a list of dictionaries. """          
        query = """ SELECT id, ui_alias FROM %s ORDER BY cores """ % self.table
        return self.conn.rs(query, how=self.conn.RETURN_DICT)

    def credentials(self, server_id):
        """ Retrieve ssh login information for a server. """ 
        query = """ 
        SELECT host, username, private_key 
        FROM %s 
        WHERE id = %s
        """ % (self.table, server_id)
        return self.conn.row(query, how=self.conn.RETURN_DICT)
    
    def path(self, server_id):
        """ Retrieve the root path where the OPALgui stores its data. """
        query = """
        SELECT path
        FROM %s
        WHERE id = %s
        """ % (self.table, server_id)
        return self.conn.field(query)
    
    def connection(self, server_id):
        """ Open a ssh connection and return the SSHConnection instance. """
        credentials = self.credentials(server_id)
        ssh = SSHConnection(credentials['host'], credentials['username'], '~/.ssh/' + credentials['private_key'])
        return ssh
        
