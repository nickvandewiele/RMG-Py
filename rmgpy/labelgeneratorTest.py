import unittest

from rmgpy.species import Species
from rmgpy.molecule import Molecule

from .labelgenerator import *

class TestIs_chemkin_compatible(unittest.TestCase):

	def testBracketForbidden(self):
		label = '[CH3]'
		self.assertFalse(is_chemkin_compatible(label))

	def testAllowed(self):
		label = 'CH3'
		self.assertTrue(is_chemkin_compatible(label))

	def testTooLongForbidden(self):
		label = 'aaaaaaaaaaaaaaaaaaaaaaa'
		self.assertFalse(is_chemkin_compatible(label))

	def testEqualsSignForbidden(self):
		label = 'C=C'
		self.assertFalse(is_chemkin_compatible(label))

	def testHashSignForbidden(self):
		label = 'C#C'
		self.assertFalse(is_chemkin_compatible(label))


class TestGenerateSpeciesIdentifier(unittest.TestCase):
	
	def testSpecies_use_SMILES(self):
		smi = 'C'
		spc = Species(molecule=[Molecule().fromSMILES(smi)])
		kwargs = {'species':spc}
		label = generateSpeciesIdentifier(**kwargs)
		self.assertEquals(label, 'C')

	def test_species_use_formula(self):
		smi = '[CH3]'
		spc = Species(molecule=[Molecule().fromSMILES(smi)])
		
		kwargs = {
			'species':spc,
			'index':100,
		}

		label = generateSpeciesIdentifier(**kwargs)
		exp = 'CH3(100)'
		self.assertEquals(label, exp)

	def test_inchispecies_use_formula_H(self):
		inchi = 'InChI=1S/H/u1'
		kwargs = {'id':inchi}
		label = generateSpeciesIdentifier(**kwargs)
		self.assertEquals(label, 'H')

	def test_inchispecies_use_formula(self):
		inchi = 'InChI=1S/CH3/h1H3/u1'
		kwargs = {
			'id':inchi,
			'index':100,
		}
		
		label = generateSpeciesIdentifier(**kwargs)
		exp = 'CH3(100)'
		self.assertEquals(label, exp)

	def testReactive(self):
		inchi = 'InChI=1S/CH4/h1H4'
		reactive = True
		kwargs = {
			'id':inchi,
			'reactive':reactive,
		}
		label = generateSpeciesIdentifier(**kwargs)
		self.assertIsNotNone(label)

	def testUserDefined(self):
		name = 'methane'
		smi = 'C'
		spc = Species(molecule=[Molecule().fromSMILES(smi)])
		kwargs = {
			'species':spc,
			'user': name,
			}
		label = generateSpeciesIdentifier(**kwargs)
		self.assertEquals(label, name)

	def testIndex(self):
		smi = 'C'
		spc = Species(molecule=[Molecule().fromSMILES(smi)])
		index = 100
		kwargs = {
			'species':spc,
			'index': index,
			}
		label = generateSpeciesIdentifier(**kwargs)
		self.assertTrue(str(index) in label, label)

if __name__ == '__main__':
	unittest.main()