from copy import deepcopy
import os
import unittest 
import logging

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.kinetics import KineticsDatabase, TemplateReaction
from rmgpy.molecule import Molecule
from rmgpy.species import Species, TransitionState
from rmgpy.reaction import Reaction
from rmgpy.statmech.conformer import Conformer
from rmgpy.kinetics import Arrhenius
from rmgpy.statmech.translation import Translation, IdealGasTranslation
from rmgpy.statmech.rotation import Rotation, LinearRotor, NonlinearRotor, KRotor, SphericalTopRotor
from rmgpy.statmech.vibration import Vibration, HarmonicOscillator
from rmgpy.statmech.torsion import Torsion, HinderedRotor
from rmgpy.statmech.conformer import Conformer
from rmgpy.kinetics import Arrhenius
from rmgpy.thermo import Wilhoit
import rmgpy.constants as constants
###################################################


class TestReaction(unittest.TestCase):
    
    def setUp(self):
        path = os.path.join(settings['database.directory'],'kinetics')
       
        
        
    def makeEthylene(self):
        return Species(
            label = 'C2H4',
            conformer = Conformer(
                E0 = (44.7127, 'kJ/mol'),
                modes = [
                    IdealGasTranslation(
                        mass = (28.0313, 'amu'),
                    ),
                    NonlinearRotor(
                        inertia = (
                            [3.41526, 16.6498, 20.065],
                            'amu*angstrom^2',
                        ),
                        symmetry = 4,
                    ),
                    HarmonicOscillator(
                        frequencies = (
                            [828.397, 970.652, 977.223, 1052.93, 1233.55, 1367.56, 1465.09, 1672.25, 3098.46, 3111.7, 3165.79, 3193.54],
                            'cm^-1',
                        ),
                    ),
                ],
                spinMultiplicity = 1,
                opticalIsomers = 1,
            ),
            molecule=[Molecule().fromSMILES('C=C')]
        )
        
        
    def makeHydrogen(self):
        return Species(          
            label = 'H',
            conformer = Conformer(
                E0 = (211.794, 'kJ/mol'),
                modes = [
                    IdealGasTranslation(
                        mass = (1.00783, 'amu'),
                    ),
                ],
                spinMultiplicity = 2,
                opticalIsomers = 1,
            ),
            molecule=[Molecule().fromSMILES('[H]')]
        )
        
        
    def makeEthyl(self):
        return Species(
            label = 'C2H5',
            conformer = Conformer(
                E0 = (111.603, 'kJ/mol'),
                modes = [
                    IdealGasTranslation(
                        mass = (29.0391, 'amu'),
                    ),
                    NonlinearRotor(
                        inertia = (
                            [4.8709, 22.2353, 23.9925],
                            'amu*angstrom^2',
                        ),
                        symmetry = 1,
                    ),
                    HarmonicOscillator(
                        frequencies = (
                            [482.224, 791.876, 974.355, 1051.48, 1183.21, 1361.36, 1448.65, 1455.07, 1465.48, 2688.22, 2954.51, 3033.39, 3101.54, 3204.73],
                            'cm^-1',
                        ),
                    ),
                    HinderedRotor(
                        inertia = (1.11481, 'amu*angstrom^2'),
                        symmetry = 6,
                        barrier = (0.244029, 'kJ/mol'),
                        semiclassical = None,
                    ),
                ],
                spinMultiplicity = 2,
                opticalIsomers = 1,
            ),
            molecule=[Molecule().fromSMILES('C[CH2]')]
        )
        
        
    def makeTS(self):
        return TransitionState(
            label = 'TS',
            conformer = Conformer(
                E0 = (266.694, 'kJ/mol'),
                modes = [
                    IdealGasTranslation(
                        mass = (29.0391, 'amu'),
                    ),
                    NonlinearRotor(
                        inertia = (
                            [6.78512, 22.1437, 22.2114],
                            'amu*angstrom^2',
                        ),
                        symmetry = 1,
                    ),
                    HarmonicOscillator(
                        frequencies = (
                            [412.75, 415.206, 821.495, 924.44, 982.714, 1024.16, 1224.21, 1326.36, 1455.06, 1600.35, 3101.46, 3110.55, 3175.34, 3201.88],
                            'cm^-1',
                        ),
                    ),
                ],
                spinMultiplicity = 2,
                opticalIsomers = 1,
            ),
            frequency = (-750.232, 'cm^-1'),
        )
    
    
    def makeKinetics(self):
        return Arrhenius(
                A = (501366000.0, 'cm^3/(mol*s)'),
                n = 1.637,
                Ea = (4.32508, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (2500, 'K'),
            )
        
    
    def makeLibraryReaction(self):
        '''C2H4 + H <=> C2H5'''
        path = os.path.join(settings['database.directory'],'kinetics','libraries', 'combustion_core')
        
        db = KineticsDatabase()
        db.libraryOrder.append(('version5','Reaction Library'))
        db.loadLibraries(path, libraries=['version5'])
        lib = db.libraries['version5']
        
        reactants = [Molecule().fromSMILES('C=C'), Molecule().fromAdjacencyList("""
        1 H u1 p0 c0 
            """
        )]
        products = [Molecule().fromSMILES('C[CH2]')]
        
        return db.generateReactionsFromLibrary(reactants, products, lib)[0]

    
    def makeTemplateReaction(self):
        '''C2H4 + H <=> C2H5'''
        path = os.path.join(settings['database.directory'],'kinetics','families')
        
        db = KineticsDatabase()
        db.loadFamilies(path, families=['R_Addition_MultipleBond'])
        
        #load forbidden structures
        database = RMGDatabase()
        database.loadForbiddenStructures(os.path.join(settings['database.directory'],'forbiddenStructures.py'))
        
        reactants = [Molecule().fromSMILES('C=C'), Molecule().fromAdjacencyList("""
        1 H u1 p0 c0 
            """
        )]
        products = [Molecule().fromSMILES('C[CH2]')]

        return db.generateReactionsFromFamilies(reactants, products, only_families='R_Addition_MultipleBond')[0]
    
    def makeGenericReaction(self):
        '''C2H4 + H <=> C2H5'''
        self.ethylene = self.makeEthylene()
        self.ethyl = self.makeEthyl()
        self.h = self.makeHydrogen()
        self.TS = self.makeTS()
        self.kinetics = self.makeKinetics()

        return Reaction(
            reactants = [self.h, self.ethylene],
            products = [self.ethyl], 
            kinetics = self.kinetics,
            transitionState = self.TS,
        )
        
    def testEquality(self):
        rxn_gen = self.makeGenericReaction()
        rxn_lib = self.makeLibraryReaction()
        rxn_fam = self.makeTemplateReaction()
        
        self.assertEqual(rxn_fam, rxn_gen)
        self.assertEqual(rxn_fam, rxn_lib)
        self.assertEqual(rxn_lib, rxn_gen)
        
        #check for symmetry
        self.assertEqual(rxn_fam == rxn_gen, rxn_gen == rxn_fam)
        self.assertEqual(rxn_fam == rxn_lib, rxn_lib == rxn_fam)
        self.assertEqual(rxn_lib == rxn_gen, rxn_gen == rxn_lib)
        
        from sets import Set
        s = Set([rxn_fam, rxn_gen, rxn_lib])#
        self.assertEqual(len(s), 1)
    
    def testDuplicateReactionsEquality(self):
        rxn_gen = self.makeGenericReaction()
        rxn_fam = self.makeTemplateReaction()
        rxn_lib = self.makeLibraryReaction()
        rxn_lib.duplicate = True
        
        from copy import deepcopy
        rxn_lib2 = deepcopy(rxn_lib)
        self.assertNotEqual(rxn_lib, rxn_lib2)
        
        from sets import Set
        s = Set([rxn_lib, rxn_lib2])
        self.assertEqual(len(s), 2)
        
        #although the library reaction has the duplicate attribute, it is considered equal to the generic reaction
        self.assertEqual(rxn_lib, rxn_gen)
        self.assertEqual(rxn_lib, rxn_fam)
        
        rxn_fam.duplicate = True
        rxn_gen.duplicate = True
        '''
        Even when generic or template reactions have a duplicate attribute, the fact that they originate from different
        sources makes them identical.
        '''
        
        self.assertEqual(rxn_lib, rxn_gen)
        self.assertEqual(rxn_lib, rxn_fam)
    
    def testNNH_GRIMech3N(self):
        '''Test that two duplicate reactions of NNH <=> H + N2 are correctly read.'''
        path = os.path.join(settings['database.directory'],'kinetics','libraries')
        label = 'GRI-Mech3.0-N'
        db = KineticsDatabase()
        db.libraryOrder.append((label,'Reaction Library'))
        db.loadLibraries(path, libraries=[label])
        lib = db.libraries[label]
        
        reactants = [ Molecule().fromAdjacencyList("""
        multiplicity 2
        1 N u0 p1 c0 {2,D} {3,S}
        2 N u1 p1 c0 {1,D}
        3 H u0 p0 c0 {1,S}
            """
        )]
        products = [Molecule().fromSMILES('N#N'), Molecule().fromAdjacencyList("""
        1 H u1 p0 c0 
            """
        )]
        
        rxns = db.generateReactionsFromLibrary(reactants, products, lib)
        self.assertEqual(len(rxns), 2)
        logging.info(rxns)
        self.assertTrue(rxns[0].duplicate)
        self.assertTrue(rxns[1].duplicate)
        self.assertFalse(rxns[0] == rxns[1])
        
        