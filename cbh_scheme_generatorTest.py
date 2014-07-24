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

class Abstract_CBH(unittest.TestCase):
                
    def check(self, d, rxn):
        shared_items = set(d.items()) & set(rxn.coefficients.items())
        self.assertEquals(len(shared_items), max(len(d), len(rxn.coefficients)))
    
    def runAbstract(self, cbh, s):
        d = read_reaction_string(s)
        cbh.run()
        rxn =  cbh.error_reaction
        self.check(d, rxn)
        
class TestCBH0(Abstract_CBH):
    def __init__(self, *args, **kwargs):
            super(self.__class__, self).__init__(*args, **kwargs)
 
    def runCBH(self,smi, s):
        spc = gen.makeSpeciesFromSMILES(smi)
        cbh = gen.CBH0Reaction(spc=spc)
        self.runAbstract(cbh, s)
        
    def testCPD(self):
        smi = 'C1=CCC=C1'
        s = '{1}C1=CCC=C1 + {7}[H][H] <=> {5}C'
        self.runCBH(smi, s)
        
    def testOxazoline(self):
        smi = 'C1CN=CO1'
        s = '{1}C1CN=CO1 + {6}[H][H] <=> {3}C + {1}N + {1}O'
        self.runCBH(smi,s)

    def testNHeptanethiol(self):
        smi = 'C(CCCC)CCS'
        s = '{1}C(CCCC)CCS + {7}[H][H] <=> {7}C + {1}S'
        self.runCBH(smi, s)
        
    def testCyclohexanone(self):
        smi = 'C1CCC(=O)CC1'
        s = '{1}C1CCC(=O)CC1 + {8}[H][H] <=> {6}C + {1}O'
        self.runCBH(smi, s)
        
class TestCBH1(Abstract_CBH):
            
    def runCBH(self,smi, s):
        spc = gen.makeSpeciesFromSMILES(smi)
        cbh = gen.CBH1Reaction(spc=spc)
        self.runAbstract(cbh, s)
    
    def testCPD(self):
        smi = 'C1C=CC=C1'
        s = '{1}C1C=CC=C1 + {5}C <=> {3}CC + {2}C=C'
        self.runCBH(smi,s)

    def testOxazoline(self):
        smi = 'C1CN=CO1'
        s = '{1}C1CN=CO1 + {3}C + {1}N + {1}O <=> {1}CC + {2}CO + {1}CN + {1}C=N'
        self.runCBH(smi, s)

    def testNHeptanethiol(self):
        smi = 'C(CCCC)CCS'
        s = '{1}C(CCCC)CCS + {6}C <=> {6}CC + {1}CS'
        self.runCBH(smi, s)
        
    def testCyclohexanone(self):
        smi = 'C1CCC(=O)CC1'
        s = '{1}C1CCC(=O)CC1 + {7}C <=> {6}CC + {1}C=O'
        self.runCBH(smi, s)
        
    def testPropylPent4Enoate(self):
        smi = 'C(COC(=O)CCC=C)C'
        s = '{1}C(COC(=O)CCC=C)C + {7}C + {1}O <=> {5}CC + {1}C=C + {1}C=O + {2}CO'
        self.runCBH(smi, s)

class TestCBH2(Abstract_CBH):
            
    def runCBH(self,smi, s):
        spc = gen.makeSpeciesFromSMILES(smi)
        cbh = gen.CBH2Reaction(spc=spc)
        self.runAbstract(cbh, s)
    
    def testCPD(self):
        smi = 'C1C=CC=C1'
        s = '{1}C1C=CC=C1 + {3}CC + {2}C=C <=> {1}CCC + {4}C=CC'
        self.runCBH(smi,s)

    def testOxazoline(self):
        smi = 'C1CN=CO1'
        s = '{1}C1CN=CO1 + {1}CC + {2}CO + {1}CN + {1}C=N <=> {1}C(C)N + {1}CCO + {1}C(=N)O + {1}CN=C + {1}COC'
        self.runCBH(smi, s)

    def testNHeptanethiol(self):
        smi = 'C(CCCC)CCS'
        s = '{1}C(CCCC)CCS + {5}CC <=> {5}CCC + {1}CCS'
        self.runCBH(smi, s)
        
    def testCyclohexanone(self):
        smi = 'C1CCC(=O)CC1'
        s = '{1}C1CCC(=O)CC1 + {6}CC <=> {5}CCC + {1}CC(C)=O'
        self.runCBH(smi, s)
        
    def testPropylPent4Enoate(self):
        smi = 'C(COC(=O)CCC=C)C'
        s = '{1}C(COC(=O)CCC=C)C + {4}CC + {2}CO <=> {3}CCC + {1}C=CC + {1}CC(=O)O + {1}COC + {1}CCO'
        self.runCBH(smi, s)
 
class TestCBH3(Abstract_CBH):
            
    def runCBH(self,smi, s):
        spc = gen.makeSpeciesFromSMILES(smi)
        cbh = gen.CBH3Reaction(spc=spc)
        self.runAbstract(cbh, s)
    
    def testCPD(self):
        smi = 'C1C=CC=C1'
        s = '{1}C1C=CC=C1 + {1}CCC + {4}C=CC <=> {2}C=CCC + {1}C=CC=C + {2}CC=CC'
        self.runCBH(smi,s)

    def testOxazoline(self):
        smi = 'C1CN=CO1'
        s = '{1}C1CN=CO1 + {1}C(C)N + {1}CCO + {1}C(=N)O + {1}CN=C + {1}COC <=> {1}C(CO)N + {1}CCOC + {1}CN=CO + {1}COC=N + {1}C(C)N=C'
        self.runCBH(smi, s)

    def testNHeptanethiol(self):
        smi = 'C(CCCC)CCS'
        s = '{1}C(CCCC)CCS + {4}CCC <=> {4}CCCC + {1}CCCS'
        self.runCBH(smi, s)
        
    def testCyclohexanone(self):
        smi = 'C1CCC(=O)CC1'
        s = '{1}C1CCC(=O)CC1 + {5}CCC + {1}CC(C)=O <=> {4}CCCC + {2}CCC(C)=O'
        self.runCBH(smi, s)
        
    def testPropylPent4Enoate(self):
        smi = 'C(COC(=O)CCC=C)C'
        s = '{1}C(COC(=O)CCC=C)C + {2}CCC + {1}COC + {1}CCO + {1}CC(=O)O <=> {1}C=CCC + {1}CCCC + {1}CCC(=O)O + {1}COC(C)=O + {1}CCOC + {1}CCCO'
        self.runCBH(smi, s)
                
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()