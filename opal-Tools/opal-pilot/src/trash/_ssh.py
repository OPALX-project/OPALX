import commands
import string
import os
import re

class ClientSSH():
    """
    A very simple and limited SSH client.
    
    Simple SSH client that operates directly on the shell. Not compatible to 
    Windows or Python 3.x. Replace by something decent later.
    """

    def __init__(self, privkey, user, host):
        self.filters = ["(module: command not found)$"]
        self.__privkey = "$HOME/.ssh/%s" % privkey
        self.__user = user
        self.__host = host

    def isdir(self, dir):
        if self.getline('if [ -d "%s" ]; then echo "true"; fi' % (dir)) == 'true':
            return True
        return False

    def isfile(self, file):
        if self.getline('if [ -f "%s" ]; then echo "true"; fi' % (file)) == 'true':
            return True
        return False

    def execute(self, command):
        self.getoutput(command)

    def getline(self, command):
        output = self.getoutput(command)
        lines = output.splitlines()

        if len(lines) >= 1:
            return lines[0]
        else:
            return None
        
    def getoutput(self, command):
        """ commands is deprecated with Python 3.0
        see: http://docs.python.org/library/commands.html
        """
        output = commands.getoutput("ssh -i %s %s@%s '%s'" % (self.__privkey, self.__user, self.__host, command))
        lines = output.splitlines()
        linesnew = []
        for line in lines:
            for filter in self.filters:
                if re.search(filter, line):
                    line = False
                    break
            if line:
                linesnew.append(line)

        return string.join(linesnew, os.linesep)
    
    def mkdir(self, directory):
        self.execute("mkdir %s" % directory)
    
    def link(self, target, directory):
        self.execute("ln -s %s %s" % (target, directory))

    def scp(self, file_src, file_target):
        """ commands is deprecated with Python 3.0
        see: http://docs.python.org/library/commands.html
        """
        commands.getoutput("scp -i %s %s %s@%s:%s" % (self.__privkey, file_src, self.__user, self.__host, file_target))