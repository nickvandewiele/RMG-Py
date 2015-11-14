import os
import sys
import unittest

from .reduction import *
from .input import load


class MockMolecule(object):
    """docstring for MockMolecule"""
    def __init__(self, label):
        super(MockMolecule, self).__init__()
        self.label = label
        
class ReductionReactionTest(unittest.TestCase):

    def setUp(self):
        from rmgpy.reaction import Reaction
        from .model import ReductionReaction

        mol1 = MockMolecule(label='mol1')
        mol2 = MockMolecule(label='mol2')
        mol3 = MockMolecule(label='mol3')
        mol4 = MockMolecule(label='mol4')
        
        self.rxn = Reaction(reactants=[mol1, mol2], products=[mol3, mol4])
        
        self.rrxn = ReductionReaction(self.rxn)


    def tearDown(self):
        del self.rrxn


    def test_constructor(self):
        rrxn = self.rrxn
        rxn = self.rxn

        self.assertIsNotNone(rrxn)

        # attributes
        self.assertIsNotNone(rrxn.reactants, rxn.reactants)
        self.assertIs(rrxn.products, rxn.products)
        self.assertIs(rrxn.rmg_reaction, rxn)
        self.assertIsNotNone(rrxn.stoichio)
        self.assertIsNone(rrxn.kf)
        self.assertIsNone(rrxn.kb)


        # stoichio
        for k,d in self.rrxn.stoichio.iteritems():
            for k,v in d.iteritems():
                self.assertEquals(v, 1)



    def test_reduce(self):
        import pickle
        reaction = pickle.loads(pickle.dumps(self.rrxn))

class OptimizeTest(unittest.TestCase):

    #MINIMAL
    working_dir = 'rmgpy/reduction/test_data/minimal/'
    inputFile = working_dir + 'input.py'
    reductionFile = working_dir + 'reduction_input.py'
    chemkinFile = working_dir + 'chemkin/chem.inp'
    spc_dict = working_dir + 'chemkin/species_dictionary.txt'


    @classmethod
    def setUpClass(cls):
        super(OptimizeTest, cls).setUpClass()
        rmg, target_label, error = load(cls.inputFile, cls.reductionFile, cls.chemkinFile, cls.spc_dict)
        cls.rmg = rmg
        cls.target_label = target_label
        cls.error = error

        reactionModel = rmg.reactionModel
        initialize(rmg.outputDirectory, reactionModel.core.reactions)
    

    def test_compute_conversion(self):
        rmg = OptimizeTest.rmg
        target_label = OptimizeTest.target_label
        reactionModel = rmg.reactionModel

        atol, rtol = rmg.absoluteTolerance, rmg.relativeTolerance
        index = 0
        reactionSystem = rmg.reactionSystems[index]

        conv = compute_conversion(target_label, reactionModel, reactionSystem, index,\
         rmg.absoluteTolerance, rmg.relativeTolerance)
        self.assertIsNotNone(conv)


    def test_reduce_compute(self):
        rmg = OptimizeTest.rmg
        target_label = OptimizeTest.target_label
        reactionModel = rmg.reactionModel


        atol, rtol = rmg.absoluteTolerance, rmg.relativeTolerance
        index = 0
        reactionSystem = rmg.reactionSystems[index]

        orig_conv = compute_conversion(target_label, reactionModel, reactionSystem, index,\
         rmg.absoluteTolerance, rmg.relativeTolerance)

        tols = [0.7, 1e-3, 1e-6]
        for tol in tols:
            conv, important_rxns = reduce_compute(tol, target_label, reactionModel, rmg, index)
            self.assertIsNotNone(conv)

    def test_optimize(self):
        rmg = OptimizeTest.rmg
        target_label = OptimizeTest.target_label
        error = OptimizeTest.error
        reactionModel = rmg.reactionModel


        atol, rtol = rmg.absoluteTolerance, rmg.relativeTolerance
        index = 0
        reactionSystem = rmg.reactionSystems[index]
        
        #compute original target conversion
        Xorig = compute_conversion(target_label, reactionModel, reactionSystem, index,\
         rmg.absoluteTolerance, rmg.relativeTolerance)

        # optimize reduction tolerance
        tol, important_rxns = optimize_tolerance(target_label, reactionModel, rmg, index, error, Xorig)
        self.assertAlmostEqual(1e-04, tol)


if __name__ == '__main__':
    unittest.main()

if __name__ == '__main__':
    unittest.main()