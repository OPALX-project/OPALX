"""
This package contains the data models of the application. 

All the business logic of the application happens here.
Database manipulation should only happen in this package.

Conventions for classes in this package:
* Naming
-- First letter always upper case (Yes: Simulation, No: simulation) 
-- Camel case, not underscores (Yes: StatisticUnit, No: Statistic_Unit, Statistic_unit)
-- Keep a strong relation to the corresponding database table.
   Examples (database table -> model): 
   simulation                  -> Simulation
   accelerator_file            -> AcceleratorFile
   simulation_statistic_column -> SimulationStatisticColumn   
* Modules
-- The first part of the model name relates in which module it is found
   Examples (model -> module):
   Simulation                -> model.simulation
   Accelerator               -> model.accelerator
   AcceleratorFile           -> model.accelerator
   SimulationStatisticColumn -> model.simulation (not model.statistic)

"""

import string
import re

from util.database import DatabaseConnection
from util.design_pattern import Singleton

class Model():
    """Provides the minimal requirements for every model.
    
    A model does at least need a database connection (self.conn) and its
    database table name (self.table). The database table name is retrieved
    from the class name for convenience. Please follow the conventions
    mentioned in the top of this document.
        
    """
    
    def __init__(self):
        self.__re_cap = re.compile('([a-z0-9]+)([A-Z]+)')
        self.conn = DatabaseConnection()
        self.table = self.__classname2tablename()
              
    def __classname2tablename(self):
        classname = str(self.__class__)
        classname = string.split(classname, '.')[-1]
        classname = self.__re_cap.sub(r'\1_\2', classname)
        classname = classname.lower()
        return classname
        
        
class CachingModel(Model):
    """Extend the base class by adding caching facilities (self.cache)."""
    
    def __init__(self):
        Model.__init__(self)
        self.cache = Cache()
        
        
class Cache():
    """Provide a simple caching mechanism."""
    
    def __init__(self):
        self.__cache = {}
        
    def find(self, key):
        if key in self.__cache.keys():
            return self.__cache[key]
        return None
    
    def put(self, key, value):
        self.__cache[key] = value  
        