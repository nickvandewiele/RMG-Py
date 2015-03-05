from copy import deepcopy
import os
import unittest 

from rmgpy import settings
from rmgpy.data.kinetics import KineticsDatabase
###################################################

class TestKineticsFamily(unittest.TestCase):
    
    def setUp(self):
        path = os.path.join(settings['database.directory'],'kinetics','families')
        database = KineticsDatabase()
        database.loadFamilies(path, families=['H_Abstraction'])
        self.f1 = database.families['H_Abstraction']
        
    def testIdenticalKineticsFamilies(self):
        """
        Test that 2 identical, but distinct KineticsFamily objects are
        not added twice to a dictionary
        """
        f1c = deepcopy(self.f1)
        
        d = {}
        d[self.f1] = 'foo'
        d[deepcopy(self.f1)] = 'foo'
        self.assertEquals(1, len(d.keys()))