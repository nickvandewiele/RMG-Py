import os
import sys
import unittest

from .input import load
from .reduction import initialize, compute_conversion

from .optimization import *

class OptimizeTest(unittest.TestCase):

    #MINIMAL
    wd = os.path.join('rmgpy/reduction/test_data/minimal/')
    inputFile = os.path.join(wd, 'input.py')
    reductionFile = os.path.join(wd, 'reduction_input.py')
    chemkinFile = os.path.join(wd, 'chemkin','chem.inp')
    spc_dict = os.path.join(wd, 'chemkin','species_dictionary.txt')

    @classmethod
    def setUpClass(cls):
        super(OptimizeTest, cls).setUpClass()
        rmg, target_label, error = load(cls.inputFile, cls.reductionFile, cls.chemkinFile, cls.spc_dict)
        cls.rmg = rmg
        cls.target_label = target_label
        cls.error = error

        reactionModel = rmg.reactionModel
        initialize(rmg.outputDirectory, reactionModel.core.reactions)
    
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
        tol, important_rxns = optimize(target_label, reactionModel, rmg, index, error, Xorig)
        self.assertAlmostEqual(1e-04, tol)


if __name__ == '__main__':
    unittest.main()
