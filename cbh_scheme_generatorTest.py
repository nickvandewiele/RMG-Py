'''
Created on Jul 21, 2014

@author: nickvandewiele
'''
import unittest
import cbh_scheme_generator as gen

class TestCBH0(unittest.TestCase):

    def runCBH0(self,smi):
        spc = gen.makeSpeciesFromSMILES(smi)
        cbh0 = gen.CBH0Reaction(spc=spc)
        cbh0.run()
        rxn =  cbh0.error_reaction
        return rxn
    
    def testCPD(self):
        smi = 'C1C=CC=C1'
        rxn = self.runCBH0(smi)
        self.assertEquals(rxn.coefficients['[H][H]'], 7)
        self.assertEquals(rxn.coefficients['C'], 5)

    def testOxazoline(self):
        smi = 'C1CN=CO1'
        rxn = self.runCBH0(smi)
        self.assertEquals(rxn.coefficients['[H][H]'], 6)
        self.assertEquals(rxn.coefficients['C'], 3)
        self.assertEquals(rxn.coefficients['N'], 1)
        self.assertEquals(rxn.coefficients['O'], 1)

    def testNHeptanethiol(self):
        smi = 'C(CCCC)CCS'
        rxn = self.runCBH0(smi)
        self.assertEquals(rxn.coefficients['[H][H]'], 7)
        self.assertEquals(rxn.coefficients['C'], 7)
        self.assertEquals(rxn.coefficients['S'], 1)
    def testCyclohexanone(self):
        smi = 'C1CCC(=O)CC1'
        rxn = self.runCBH0(smi)
        self.assertEquals(rxn.coefficients['[H][H]'], 8)
        self.assertEquals(rxn.coefficients['C'], 6)
        self.assertEquals(rxn.coefficients['O'], 1)
        
class TestCBH1(unittest.TestCase):

    def runCBH1(self,smi):
        spc = gen.makeSpeciesFromSMILES(smi)
        cbh1 = gen.CBH1Reaction(spc=spc)
        cbh1.run()
        rxn =  cbh1.error_reaction
        return rxn
    
    def testCPD(self):
        smi = 'C1C=CC=C1'
        rxn = self.runCBH1(smi)
        self.assertEquals(rxn.coefficients['C'], 5)
        
        self.assertEquals(rxn.coefficients['CC'], 3)
        self.assertEquals(rxn.coefficients['C=C'], 2)

    def testOxazoline(self):
        smi = 'C1CN=CO1'
        rxn = self.runCBH1(smi)
        self.assertEquals(rxn.coefficients['C'], 3)
        self.assertEquals(rxn.coefficients['N'], 1)
        self.assertEquals(rxn.coefficients['O'], 1)
        
        self.assertEquals(rxn.coefficients['CC'], 1)
        self.assertEquals(rxn.coefficients['CO'], 2)
        self.assertEquals(rxn.coefficients['CN'], 1)
        self.assertEquals(rxn.coefficients['C=N'], 1)

    def testNHeptanethiol(self):
        smi = 'C(CCCC)CCS'
        rxn = self.runCBH1(smi)
        self.assertEquals(rxn.coefficients['C'], 6)
        
        self.assertEquals(rxn.coefficients['CC'], 6)
        self.assertEquals(rxn.coefficients['CS'], 1)
        
    def testCyclohexanone(self):
        smi = 'C1CCC(=O)CC1'
        rxn = self.runCBH1(smi)
    
        self.assertEquals(rxn.coefficients['C'], 7)
        
        self.assertEquals(rxn.coefficients['CC'], 6)
        self.assertEquals(rxn.coefficients['C=O'], 1)    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()