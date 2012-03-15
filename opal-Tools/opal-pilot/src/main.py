"""
Provide a GUI for OPAL with provenance in mind.

This GUI provides a easy way to get OPAL simulations started that are based
on nearly static settings. It is based on the MVC pattern.

MVC = Model-View-Controller
Models contain the business logic for the application and provide way to access
and manipulate data.
Controllers accept requests and work with models to prepare data for the view.
Views take data from the model and controller and present it in a way the
user/client understands.

@author: David Uhlig <david.uhlig@googlemail.com>
@organization: Paul Scherrer Institut
@version: 0.8.0
  
"""
import sys

from PyQt4.QtGui import QApplication
from PyQt4.QtCore import QCoreApplication

from controller.main import MainWindowController
       
def main():
    QCoreApplication.setOrganizationName('Paul Scherrer Institut')
    QCoreApplication.setOrganizationDomain('psi.ch')
    QCoreApplication.setApplicationName('OPAL Pilot')
    
    app = QApplication(sys.argv)
    window = MainWindowController()
    sys.exit(app.exec_())
    
if __name__ == "__main__":
    main()
