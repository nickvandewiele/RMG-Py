from copy import deepcopy
import os
import unittest 

from rmgpy import settings
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

class TestTemplateReaction(unittest.TestCase):
    
    def setUp(self):
        path = os.path.join(settings['database.directory'],'kinetics','families')
        self.database = KineticsDatabase()
        self.database.loadFamilies(path, families=['H_Abstraction', '1,3_Insertion_CO2'])
        
        self.ethylene = self.makeEthylene()
        self.ethyl = self.makeEthyl()
        self.h = self.makeHydrogen()
        self.TS = self.makeTS()
        self.kinetics = self.makeKinetics()
        
        self.reaction1 = self.makeReaction1()
        self.reaction2 = self.makeReaction2()
        
        
        
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
        
        
    def makeReaction1(self):
        f1 = self.database.families['H_Abstraction']

        return TemplateReaction(
            reactants = [self.h, self.ethylene],
            products = [self.ethyl], 
            kinetics = self.kinetics,
            transitionState = self.TS,
            family = f1
        )
        
    def makeReaction1Reverse(self):
        f1 = self.database.families['H_Abstraction']

        return TemplateReaction(
            reactants = [self.ethyl], 
            products = [self.h, self.ethylene],
            kinetics = None,
            transitionState = self.TS,
            family = f1
        )
        
        
    def makeReaction2(self):
        f2 = self.database.families['1,3_Insertion_CO2']
        
        return TemplateReaction(
            reactants = [self.h, self.ethylene],
            products = [self.ethyl], 
            kinetics = self.kinetics,
            transitionState = self.TS,
            family = f2
        )
        
    def testEquality(self):
        from copy import deepcopy
        
        self.assertNotEqual(self.reaction1.family, self.reaction2.family)
        
        self.assertFalse(self.reaction1.family == self.reaction2.family)
        self.assertTrue(self.reaction1.family != self.reaction2.family)
        
        self.assertTrue(self.reaction1.family == deepcopy(self.reaction1.family))
        self.assertFalse(self.reaction1.family != deepcopy(self.reaction1.family))
        
        self.assertFalse(self.reaction1 == self.reaction2)
        self.assertTrue(self.reaction1 != self.reaction2)

        self.assertTrue(self.reaction1 == deepcopy(self.reaction1))
        self.assertFalse(self.reaction1 != deepcopy(self.reaction1))
        
    def testEqualityReverseReaction(self):
        r1 = self.makeReaction1()
        r1v = self.makeReaction1Reverse()
        
        self.assertEqual(r1, r1v)
        
        pass
    
    def testCopy(self):
        r1 = self.makeReaction1()
        r1.labeledAtoms = 'foo'
        other = r1.copy()
        self.assertEquals(r1.labeledAtoms, other.labeledAtoms)