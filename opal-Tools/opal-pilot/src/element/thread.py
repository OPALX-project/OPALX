import os
import re
import warnings
from threading import Thread

from model.server import Server
from model.simulation import Simulation
from model.statistic import Statistic
from util.xml_parser import QStatParser

class Postprocessing(Thread):
    """Postprocess one single simulation.
    
    At this point postprocessing only contains parsing
    the .stat file and write the values to the database.
    
    """
    
    def run(self):       
        server = Server()
        simulation = Simulation()
        statistic = Statistic()
        next = simulation.nextToParse()
        
        if next is None:
            return
          
        simulation_id = next['id']
        server_id = next['server_id']
        process_id = next['process_id']
        
        simulation.setState(process_id, Simulation.STATE_POSTPROCESSING)
        
        ssh = server.connection(server_id)
        path = '%s/%s/opalgui.stat' % (server.path(server_id), simulation_id)
        
        if not ssh.isfile(path):
            simulation.setState(process_id, Simulation.STATE_ABORTED)
            return
        
        localpath = 'tmp/stat.%s' % simulation_id
        if not os.path.exists(localpath):
            open(localpath, 'wb').close()        
        ssh.get(path, localpath)
        
        with open(localpath, 'rb') as statfile:
            content = statfile.read()
            statistic.parseStatFile(simulation_id, content)
            
        simulation.setState(process_id, Simulation.STATE_FINISHED)
        

class PollServer(Thread):        
    
    def run(self):
        simulation = Simulation()
        unfinished = simulation.unfinished()
        server = Server()
        
        # Get servers that have unfinished simulations.
        servers2poll = []
        for item in unfinished.values():
            servers2poll.append(item['server_id'])
        servers2poll = set(servers2poll)
        
        for id in servers2poll:
            ssh = server.connection(id)
            path = server.path(id)
            ssh.execute('cd %s && qstat -xml > qstat.xml' % path)
            file = ssh.file('%s/qstat.xml' % path, 'r')
            xml = file.read()
            parser = QStatParser(xml)
            states = parser.states()
            
            # Update simulations that were found in the qstat output. 
            for item in states:
                process_id, state = item
                # The qstat output might have process not started by the gui.
                if not process_id in unfinished.keys():
                    continue
                if state == 'r':
                    simulation.setState(process_id, simulation.STATE_RUNNING)
                    del unfinished[process_id]
                # Still queued.
                else:
                    del unfinished[process_id]
            
            # Update the rest of the simulations.
            finished = []
            for process_id in unfinished.keys():
                item = unfinished[process_id]
                if item['server_id'] != id or process_id is None:
                    continue
                simulation_id = item['id']
                output = ssh.execute("""less %(path)s/%(simulation_id)s/opalgui.o%(process_id)s | grep 'End of input stream "opalgui.in".'""" % {'path': path, 'process_id': process_id, 'simulation_id': simulation_id})
                if len(output) is 0:
                    simulation.setState(process_id, simulation.STATE_ABORTED)
                elif re.search('No such file', output[0]) is not None:
                    simulation.setState(process_id, simulation.STATE_ABORTED)
                elif re.search('End of input stream', output[0]) is not None:
                    simulation.setState(process_id, simulation.STATE_POSTPROCESSING_QUEUE)
                else:
                    warnings.warn("Don't know what to do with output: '%s' for server_id %s, process_id %s" % (output[0], id, process_id))
                finished.append(process_id)
                
            for process_id in finished:
                del unfinished[process_id]          