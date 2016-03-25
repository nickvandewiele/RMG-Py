#!/usr/bin/env python
# encoding: utf-8 -*-

"""
This module contains unit tests of the rmgpy.parallel module.
"""

import os
import sys
import unittest
import random
import itertools
import logging

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.rmg.main import RMG
from rmgpy.scoop_framework.framework import TestScoopCommon
from rmgpy.molecule.molecule import Molecule
from rmgpy.data.rmg import getDB
from rmgpy.scoop_framework.util import map_, WorkerWrapper
from rmgpy.rmg.react import generate

try:
    from scoop import futures, _control, shared
except ImportError, e:
    import logging as logging
    logging.debug("Could not properly import SCOOP.")


TESTFAMILY = 'H_Abstraction'

def load():
    tearDown()
    rmg = RMG()#for solvent
    database = RMGDatabase()
    database.loadThermo(os.path.join(settings['database.directory'], 'thermo'))
    database.loadTransport(os.path.join(settings['database.directory'], 'transport'))
    path = os.path.join(settings['database.directory'])

    # forbidden structure loading
    database.loadForbiddenStructures(os.path.join(path, 'forbiddenStructures.py'))
    # kinetics family loading
    database.loadKinetics(os.path.join(path, 'kinetics'),
                                   kineticsFamilies=[TESTFAMILY],
                                   reactionLibraries=[]
                                   )
def tearDown():
    """
    Reset the loaded database
    """
    import rmgpy.data.rmg
    rmgpy.data.rmg.database = None


def funcGenerate():
    """
    Test that reaction generation from the available families works.
    """
    load()

    molA = Molecule().fromSMILES('CCCCCCCCCCCCC1C=CC=CC=1')
    molB = Molecule().fromSMILES('[OH]')

    families = getDB('kinetics').families
    family = families[TESTFAMILY]

    template = family.forwardTemplate
    forward = True
    reactantStructures = [molA, molB]
    # Reactants stored as A + B
    mappingsA = family.matchReactantToTemplate(molA, template.reactants[0])
    mappingsA = [{0: m} for m in mappingsA]
    mappingsB = family.matchReactantToTemplate(molB, template.reactants[1])
    mappingsB = [{1: m} for m in mappingsB]

    # Iterate over each pair of matches (A, B)
    mappings = list(itertools.product(mappingsA, mappingsB))
    N = len(mappings)

    results = map_(
                WorkerWrapper(generate),
                mappings,
                [reactantStructures]*N,
                [forward]*N,
                [TESTFAMILY]*N
            )
    rxnList = []
    rxnList.extend(filter(None, list(results)))

    return len(rxnList) > 0


class ParallelReactTest(TestScoopCommon):

    def __init__(self, *args, **kwargs):
        # Parent initialization
        super(self.__class__, self).__init__(*args, **kwargs)
        
        # Only setup the scoop framework once, and not in every test method:
        super(self.__class__, self).setUp()

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "test currently only runs on linux")
    def testGenerate(self):
        """
        Test that reactions can be generated in parallel for a list of given mappings.
        """
        result = futures._startup(funcGenerate)
        self.assertEquals(result, True)

if __name__ == '__main__' and os.environ.get('IS_ORIGIN', "1") == "1":
    unittest.main()