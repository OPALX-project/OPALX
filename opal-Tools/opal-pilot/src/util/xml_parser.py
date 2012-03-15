import re
from xml.dom import minidom

class QStatParser:
    
    def __init__(self, xml_string):
        self.__dom = minidom.parseString(xml_string)
        
    def states(self):
        """ Parse dom and return list of tuples containing process id and state """
        states = []
        jobs = self.__dom.getElementsByTagName('job_list')
        for job in jobs:
            job_num = self.__getText(job.getElementsByTagName('JB_job_number'))
            state = self.__getText(job.getElementsByTagName('state'))
            states.append((job_num, state))
        return states
            
    def __getText(self, nodelist):
        rc = []
        for node in nodelist:
            pattern = "<%s>(.*)</%s>" % (node.tagName, node.tagName)
            data = re.search(pattern, node.toxml()).group(1)
            rc.append(data)
        return ''.join(rc)
        
        