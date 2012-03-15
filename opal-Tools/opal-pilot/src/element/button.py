from PyQt4.QtGui import QToolButton, QTreeWidget

import controller

class AcceleratorButton(QToolButton):
    
    """
    Enable dropping for QToolButton and special accelerator related behaviour.
    """
    
    EVENT_NEW = 1
    EVENT_PLOT = 2
    
    DROPS_ACCEPT_NONE = 1       # do not accept drops
    DROPS_ACCEPT = 2            # accept dropping of one single element
    DROPS_ACCEPT_MULTIPLE = 3   # accept dropping of one or more elements
    
    def __init__(self, parent, event_type=EVENT_NEW, drops=DROPS_ACCEPT_NONE):
        if not isinstance(parent, controller.main.AcceleratorController):
            raise Exception('parent needs to be a instance of AcceleratorController')
        if event_type not in (self.EVENT_NEW, self.EVENT_PLOT, None):
            raise Exception('unknown event type %s' % event_type)
        if drops not in (self.DROPS_ACCEPT_NONE, self.DROPS_ACCEPT, self.DROPS_ACCEPT_MULTIPLE):
            raise Exception('unknown drops type')
        QToolButton.__init__(self)
        self.__parent = parent
        self.__event_type = event_type
        self.__drops_policy = drops
        self.setAcceptDrops(True)
        
    def dragEnterEvent(self, event):
        if not isinstance(event.source(), QTreeWidget):
            event.ignore()
        elif self.__isDropPolicy(self.DROPS_ACCEPT_NONE):
            event.ignore()
        elif self.__isDropPolicy(self.DROPS_ACCEPT) and self.__numDropped(event) > 1:
            event.ignore()        
        else:
            event.accept()
            
    def __isDropPolicy(self, policy):
        if self.__drops_policy == policy:
            return True
        return False
            
    def __numDropped(self, event):
        return len(event.source().selectedItems())        
        
    def dropEvent(self, event):
        if self.__event_type is self.EVENT_NEW:
            tree_widget = event.source()
            items = tree_widget.selectedItems()
            self.__parent.openInputMask(items[0])
        elif self.__event_type is self.EVENT_PLOT:
            tree_widget = event.source()
            items = tree_widget.selectedItems()
            self.__parent.openPlottingWindow(items)
