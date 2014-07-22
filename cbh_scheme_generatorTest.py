'''
Created on Jul 21, 2014

@author: nickvandewiele
'''
import unittest
import cbh_scheme_generator as gen


class TestCBH0(unittest.TestCase):

    def runCBH0(self,spc):
        cbh0 = gen.CBH0Reaction(spc=spc)
        cbh0.run()
        rxn =  cbh0.error_reaction
        return rxn
    
    def testCPD(self):
        spc = gen.makeSpecies('C1C=CC=C1')
        rxn = self.runCBH0(spc)
        self.assertEquals(rxn.coefficients['[H][H]'], 7)
        self.assertEquals(rxn.coefficients['C'], 5)

    def testOxazoline(self):
        spc = gen.makeSpecies('C1CN=CO1')
        rxn = self.runCBH0(spc)
        self.assertEquals(rxn.coefficients['[H][H]'], 6)
        self.assertEquals(rxn.coefficients['C'], 3)
        self.assertEquals(rxn.coefficients['N'], 1)
        self.assertEquals(rxn.coefficients['O'], 1)

    def testNHeptanethiol(self):
        spc = gen.makeSpecies('C(CCCC)CCS')
        rxn = self.runCBH0(spc)
        self.assertEquals(rxn.coefficients['[H][H]'], 7)
        self.assertEquals(rxn.coefficients['C'], 7)
        self.assertEquals(rxn.coefficients['S'], 1)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()