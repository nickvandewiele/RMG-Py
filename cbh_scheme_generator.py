'''
Created on Jul 21, 2014

@author: nickvandewiele
'''

import os

from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import CoreEdgeReactionModel

def makeSpeciesFromMolecule(mol):
    spc, isNew = CoreEdgeReactionModel().makeNewSpecies(mol)
    return spc

def makeSpeciesFromSMILES(smi):
    mol = Molecule().fromSMILES(smi)
    return makeSpeciesFromMolecule(mol)

class ErrorCancellingReaction(Reaction):
    '''
    Adds a coefficients dictionary as an extra attribute to
    the Reaction class that stores the number of occurrences
    of each molecule in the reaction.
    
    Assumes that all molecules in the reactants/products
    only occur ones!
    
    '''
    def __init__(self, *args, **kwargs):
            super(self.__class__, self).__init__(*args, **kwargs)
            self.coefficients = {}

class Abstract_CBH_Reaction(object):
    '''
    Superclass to all implementations of the 
    Connectivity-based Hierarchy error-canceling reaction classes.
    
    '''
    def __init__(self, spc):
        '''
        The constructor takes a Species object.
        '''
        self.spc = spc
        self.error_reaction = None    
    
    def initialize(self):
        '''
        Creates an instance of ErrorCancellingReaction with 
        empty lists for reactants and products.
        
        Adds the species attribute to the list of reactants of the
        error-cancelling reaction and updates the coefficients
        attribute of that object.
        '''
        self.error_reaction = ErrorCancellingReaction(reactants = [], products = [])
        self.error_reaction.reactants.append(self.spc)
        self.error_reaction.coefficients[self.spc.label] = 1
    
    def run(self):
        '''
        The main method that carries out 
        the different steps to create the error-cancelling reaction.
        '''
        self.initialize()
        self.populate_products()
        self.populate_reactants()
        
class CBH0Reaction(Abstract_CBH_Reaction):
    '''
    Creates rung '0' of the CBH method for the creation
    of an potential error-canceling reaction.
    '''  
    def __init__(self, spc):
        super(self.__class__, self).__init__(spc)
        
    def populate_products(self):
        '''
        Simplest, atom-centric rung where each heavy atom is extracted
        in its saturated valence state (e.g., each C as CH4).
        
        for each heavy atom, a Species will be created and added to a list, e.g.:
        [CH4, CH4, H20, NH3, H2O]
        
        Next, this list is converted into a dictionary with:
            keys = the unique species created
            values = the number of occurrences of each of the unique species
        
        This is done using the collections.Counter utility and 
        requires hashable Species objects.
        '''
        spc_list = []
        molecule = self.spc.molecule[0]
        for atom in molecule.atoms:#iterate over all atoms!
            element = atom.symbol
            if not element == 'H':#only heavy atoms
                product = makeSpeciesFromSMILES(element)
                spc_list.append(product)
        
        #populate products list of the reaction, and update coefficient for each unique product:
        from collections import Counter
        unique_spc = Counter(spc_list)
        for product, no_occurr in unique_spc.iteritems():
            self.error_reaction.products.append(product)#add to products list
            self.error_reaction.coefficients[product.label] = no_occurr #update coefficient
            
    def populate_reactants(self):
        '''
        The number of hydrogen molecules equals the number of 
        covalent bonds between heavy atoms (counting a double bond
        as two covalent bonds, etc.)
        '''
        
        #right now, i don't how to elegantly convert bond orders into ordinals.
        bond_orders = {
                       'S': 1,
                       'D': 2,
                       'T': 3
                       }
        
        molecule = self.spc.molecule[0]
        balancing_dihydrogen = 0
        
        #iterate over all unique non-hydrogen bonds of the Molecule:
        molecule.sortAtoms()#don't know if this is necessary.
        for atom1 in molecule.atoms:
            for atom2 in molecule.atoms:
                if not atom1.symbol == 'H' and not atom2.symbol == 'H':#only bonds between heavy atoms
                    if molecule.hasBond(atom1, atom2):
                        if atom1.sortingLabel < atom2.sortingLabel:
                            bond = molecule.getBond(atom1, atom2)
                            balancing_dihydrogen += bond_orders[bond.order]
    
        h2 = makeSpeciesFromSMILES('[H][H]')
        self.error_reaction.reactants.append(h2)
        self.error_reaction.coefficients[h2.label] = balancing_dihydrogen
    
    
    def populate_products_alt(self):
        '''
        An alternative, and less generic implementation of the method
        that creates products according to the CBH-O rung.
        '''
        #create an element map with keys elements in the mol and values = no. of occurrences of that element
        element_map = {}
        supported_elements = 'CHONS'
        
        for element in supported_elements:
            number = self.spc.molecule[0].getNumAtoms(element)
            if number != 0:
                element_map[element] = number 
        
        #iterate over the heavy atoms and create 'saturated' molecules:
        for element, no_occurr in element_map.iteritems():
            if not element == 'H': #heavy atoms
                #create molecule from the element symbol (eg: 'C' becomes molecule methane)
                product = makeSpeciesFromSMILES(element) 
                self.error_reaction.products.append(product)#add to products list
                self.error_reaction.coefficients[product.label] = no_occurr #update coefficient
    
    
    def populate_reactants_alt(self):
        '''
        An alternative, and less generic implementation of the method
        that creates reactants according to the CBH-O rung.
        '''
        #count no. of hydrogens in products:
        h_count_prod = 0
        for prod in self.error_reaction.products:
            h_count_prod += prod.molecule[0].getNumAtoms('H') * self.error_reaction.coefficients[prod.label]
        
        #count no. of H2 molecules required as reactant for the H-element balance of the reaction:
        balancing_dihydrogen = ( h_count_prod - self.spc.molecule[0].getNumAtoms('H') ) / 2
        h2 = makeSpecies('[H][H]')
        self.error_reaction.reactants.append(h2)
        self.error_reaction.coefficients[h2.label] = balancing_dihydrogen 
   

if __name__ == '__main__':
    spc = makeSpecies('C1C=CC=C1')
    #mol = makeSpecies('C')
    cbh0 = CBH0Reaction(spc=spc)
    cbh0.run()
    rxn =  cbh0.error_reaction
    print rxn.coefficients