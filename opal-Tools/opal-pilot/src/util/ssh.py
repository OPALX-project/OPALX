""" 
Friendly Python SSH2 interface. 

Base is taken from: 
http://commandline.org.uk/python/sftp-python-really-simple-ssh/

Extended with methods:
file
isfile
link
"""

import os
import tempfile
import string

import paramiko

class SSHConnection(object):
    """Connects and logs into the specified hostname. 
    Arguments that are not given are guessed from the environment."""
    
    __DEBUG_COMMANDS = False 

    def __init__(self, host, username = None, private_key = None, password = None, port = 22):
        self.__cache_listdir = {}
        self._sftp_live = False
        self._sftp = None
        if not username:
            username = os.environ['LOGNAME']

        # Log to a temporary file.
        templog = tempfile.mkstemp('.txt', 'ssh-')[1]
        paramiko.util.log_to_file(templog)

        # Begin the SSH transport.
        self._transport = paramiko.Transport((host, port))
        self._tranport_live = True
        # Authenticate the transport.
        if password:
            # Using Password.
            self._transport.connect(username = username, password = password)
        else:
            # Use Private Key.
            if not private_key:
                # Try to use default key.
                if os.path.exists(os.path.expanduser('~/.ssh/id_rsa')):
                    private_key = '~/.ssh/id_rsa'
                elif os.path.exists(os.path.expanduser('~/.ssh/id_dsa')):
                    private_key = '~/.ssh/id_dsa'
                else:
                    raise TypeError, "You have not specified a password or key."

            private_key_file = os.path.expanduser(private_key)
            rsa_key = paramiko.RSAKey.from_private_key_file(private_key_file)
            self._transport.connect(username = username, pkey = rsa_key)
    
    def _sftp_connect(self):
        """Establish the SFTP connection."""
        if not self._sftp_live:
            self._sftp = paramiko.SFTPClient.from_transport(self._transport)
            self._sftp_live = True

    def get(self, remotepath, localpath=None):
        """Copies a file between the remote host and the local host."""
        if not localpath:
            localpath = os.path.split(remotepath)[1]
        self._sftp_connect()
        self._sftp.get(remotepath, localpath)

    def put(self, localpath, remotepath=None):
        """Copies a file between the local host and the remote host."""
        if not remotepath:
            remotepath = os.path.split(localpath)[1]
        self._sftp_connect()
        self._sftp.put(localpath, remotepath)
        
    def file(self, filename, mode='r', bufsize=-1):
        """Open a remote file and return its file instance.
        
        Be aware that remote reading from big files is extremely
        slow. Rather copy the file to a tempfile and read it out
        locally.

        """ 
        self._sftp_connect()
        return self._sftp.file(filename, mode, bufsize)
        
    def execute(self, command):
        """Execute the given commands on a remote machine."""
        if self.__DEBUG_COMMANDS:
            print command
        channel = self._transport.open_session()
        channel.exec_command(command)
        output = channel.makefile('rb', -1).readlines()
        if output:
            return output
        else:
            return channel.makefile_stderr('rb', -1).readlines()
        
    def isfile(self, path, cached=True):
        self._sftp_connect()
        filename = os.path.basename(path)
        dir = os.path.dirname(path)
        if cached is True and dir in self.__cache_listdir:
            if filename in self.__cache_listdir[dir]:
                return True
            return False
        self.__cache_listdir[dir] = self._sftp.listdir(dir)
        if filename in self.__cache_listdir[dir]:
            return True
        return False         
    
    def link(self, target, directory):
        self.execute('ln -s %s %s' % (target, directory))

    def close(self):
        """Closes the connection and cleans up."""
        # Close SFTP Connection.
        if self._sftp_live:
            self._sftp.close()
            self._sftp_live = False
        # Close the SSH Transport.
        if self._tranport_live:
            self._transport.close()
            self._tranport_live = False

    def __del__(self):
        """Attempt to clean up if not explicitly closed."""
        self.close()
