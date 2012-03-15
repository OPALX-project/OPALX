import os
import commands

class OpalDocumentation:
    
    #FIXME: check for LATEX
    def __init(self):
        self.build()

    """
    build doc
    building docs and doxygen should not be part of the regression tests
    """
    def build():
        curdir = os.getcwd()
        os.chdir("/scratch2/amas/l_felsimsvn/src/opal/opal-Doc/doc/OPAL/user_guide")
        commands.getoutput("svn update")
        commands.getoutput("make")
        commands.getoutput("makeindex opal_user_guide")
        commands.getoutput("make")
        commands.getoutput("cp opal_user_guide.pdf /afs/psi.ch/project/amas/www/docs/opal/develop.pdf")
        os.chdir(curdir)


class OpalDoxygen:

    #FIXME: check for doxygen
    def __init__(self):
        self.build()

    """
    build doxygen
    building docs and doxygen should not be part of the regression tests
    """
    def build(self):
        curdir = os.getcwd()
        os.chdir("/scratch2/amas/l_felsimsvn/src/opal/")
        commands.getoutput("doxygen")
        commands.getoutput("cp -r doc/html /afs/psi.ch/project/amas/www/docs/opal/")
        os.chdir(curdir)
