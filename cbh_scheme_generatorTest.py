'''
Created on Jul 21, 2014

@author: nickvandewiele
'''
import unittest
import cbh_scheme_generator as gen
import re

def read_reaction_string(s):
        d = {}
        
        reaction_delimiter = '<=>'
        species_delimiter = '\+'#regex: the '+' sign
        index_delimiter = '[{}]'#regex: any of the two curly brackets
        
        sides_of_reaction = [x.strip() for x in re.split(reaction_delimiter,s)]
        
        species = []
        for side in sides_of_reaction:
            species.extend([x.strip() for x in re.split(species_delimiter, side)])
        
        for spc in species:
            index, label = filter(None, re.split(index_delimiter,spc))#ignore empty elements 
            d[label] = int(index)

        return d
    
class TestCBH0(unittest.TestCase):
    
    def check(self, d, rxn):
        for label, index in d.iteritems():
            self.assertEquals(rxn.coefficients[label], index)
            
    def runCBH0(self,smi, s):
        spc = gen.makeSpeciesFromSMILES(smi)
        cbh0 = gen.CBH0Reaction(spc=spc)
        cbh0.run()
        rxn =  cbh0.error_reaction
        d = read_reaction_string(s)   
        self.check(d, rxn)
    
    def testCPD(self):
        smi = 'C1=CCC=C1'
        s = '{1}C1=CCC=C1 + {7}[H][H] <=> {5}C'
        self.runCBH0(smi, s)
        
    def testOxazoline(self):
        smi = 'C1CN=CO1'
        s = '{1}C1CN=CO1 + {6}[H][H] <=> {3}C + {1}N + {1}O'
        self.runCBH0(smi,s)

    def testNHeptanethiol(self):
        smi = 'C(CCCC)CCS'
        s = '{1}C(CCCC)CCS + {7}[H][H] <=> {7}C + {1}S'
        self.runCBH0(smi, s)
        
    def testCyclohexanone(self):
        smi = 'C1CCC(=O)CC1'
        s = '{1}C1CCC(=O)CC1 + {8}[H][H] <=> {6}C + {1}O'
        self.runCBH0(smi, s)
        
class TestCBH1(unittest.TestCase):
    def check(self, d, rxn):
        for label, index in d.iteritems():
            self.assertEquals(rxn.coefficients[label], index)
            
    def runCBH1(self,smi, s):
        spc = gen.makeSpeciesFromSMILES(smi)
        cbh1 = gen.CBH1Reaction(spc=spc)
        cbh1.run()
        rxn =  cbh1.error_reaction
        d = read_reaction_string(s)   
        self.check(d, rxn)
    
    def testCPD(self):
        smi = 'C1C=CC=C1'
        s = '{1}C1C=CC=C1 + {5}C <=> {3}CC + {2}C=C'
        self.runCBH1(smi,s)

    def testOxazoline(self):
        smi = 'C1CN=CO1'
        s = '{1}C1CN=CO1 + {3}C + {1}N + {1}O <=> {1}CC + {2}CO + {1}CN + {1}C=N'
        self.runCBH1(smi, s)

    def testNHeptanethiol(self):
        smi = 'C(CCCC)CCS'
        s = '{1}C(CCCC)CCS + {6}C <=> {6}CC + {1}CS'
        self.runCBH1(smi, s)
        
    def testCyclohexanone(self):
        smi = 'C1CCC(=O)CC1'
        s = '{1}C1CCC(=O)CC1 + {7}C <=> {6}CC + {1}C=O'
        self.runCBH1(smi, s)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()