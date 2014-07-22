'''
Created on Jul 21, 2014

@author: nickvandewiele
'''

import os

from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import CoreEdgeReactionModel

def makeSpeciesFromMolecule(mol, label=''):
    spc, isNew = CoreEdgeReactionModel().makeNewSpecies(mol, label=label)
    return spc

def makeSpeciesFromSMILES(smi):
    mol = Molecule().fromSMILES(smi)
    return makeSpeciesFromMolecule(mol, label=smi)

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
            
            
class CBHSpeciesGenerator(object):
    def __init__(self):
        pass
    def createAdjacencyList_cbh1_product(self, atom1, atom2, bond):
        order = bond.order
        first_line = '1 '+atom1.symbol+' '+str(atom1.radicalElectrons)+' '+ str(atom1.lonePairs) +' {2,'+order+'}'
        second_line = '2 '+atom2.symbol+' '+str(atom2.radicalElectrons)+' '+ str(atom2.lonePairs) +' {1,'+order+'}'
        adjList = first_line + '\n' + second_line
        
        return adjList
    
    def createAdjacencyList_cbh2_product(self, atom, neighbors, molecule):
        lines = []
        first_line = ' '.join(['1',atom.symbol,str(atom.radicalElectrons),str(atom.lonePairs)])
        connectivity = ''
        for i, neigh in enumerate(neighbors):
            order = molecule.getBond(atom, neigh).order
            index_neigh = str(i+2)#neighbours start at index 2
            connectivity += '{'+index_neigh+','+order+'} '
            line = ' '.join([index_neigh,neigh.symbol,str(neigh.radicalElectrons),str(neigh.lonePairs),'{1,'+order+'}'])
            lines.append(line)
            
        first_line += ' '+connectivity
        lines.insert(0, first_line)
        
        return '\n'.join(lines)
    
    def create_cbh0_product(self, atom):
        return makeSpeciesFromSMILES(atom.symbol)
        
    def create_cbh1_product(self, atom1, atom2, bond):
        adjList = self.createAdjacencyList_cbh1_product(atom1, atom2, bond)#possibly a de-tour by creating adjList
        mol = Molecule().fromAdjacencyList(adjList, saturateH=True)#sature with hydrogens
        product = makeSpeciesFromMolecule(mol)
        return product 
    
    def create_cbh2_product(self, atom, neighbors, molecule):
        adjList = self.createAdjacencyList_cbh2_product(atom, neighbors, molecule)#possibly a de-tour by creating adjList
        mol = Molecule().fromAdjacencyList(adjList, saturateH=True)#sature with hydrogens
        product = makeSpeciesFromMolecule(mol)
        return product
    
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
        
    def map_species_list(self, reactants_products_list, spc_list):
        #populate reactants or products list of the reaction, and update coefficient for each unique species:
        from collections import Counter
        unique_spc = Counter(spc_list)
        for unique_species, no_occurr in unique_spc.iteritems():
            reactants_products_list.append(unique_species)#add to products list
            self.error_reaction.coefficients[unique_species.label] = no_occurr #update coefficient
    
    def exclude_hydrogens(self, atoms):
        return [atom for atom in atoms if not atom.symbol == 'H']
    
    


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
        heavy_atoms = self.exclude_hydrogens(molecule.atoms)
        for atom in heavy_atoms:#iterate over all atoms!
            product = makeSpeciesFromSMILES(atom.symbol)
            spc_list.append(product)
        
        self.map_species_list(self.error_reaction.products, spc_list)
            
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
        for atom1 in molecule.vertices:
            for atom2 in atom1.edges:
                if not atom1.symbol == 'H' and not atom2.symbol == 'H':#only bonds between heavy atoms
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
        h2 = makeSpeciesFromSMILES('[H][H]')
        self.error_reaction.reactants.append(h2)
        self.error_reaction.coefficients[h2.label] = balancing_dihydrogen 

class CBH1Reaction(Abstract_CBH_Reaction):
    '''
    Creates rung '1' of the CBH method for the creation
    of an potential error-canceling reaction.
    
    Corresponds to the Pople's isodesmic bond separation scheme.
    '''  
    def __init__(self, spc):
        super(self.__class__, self).__init__(spc)
        
        
    def populate_products(self):  
        '''
        Extracting all heavy-atom bonds in the molecule as 
        isolated valence-satisfied molecules 
        (e.g. C-C double bond as H2C=CH2).
        '''

        spc_list = []

        #iterate over all heavy-atom bonds:
        molecule = self.spc.molecule[0]
        
        #iterate over all unique non-hydrogen bonds of the Molecule:
        molecule.sortAtoms()#don't know if this is necessary.
        for atom1 in molecule.vertices:
            for atom2 in atom1.edges:
                if not atom1.symbol == 'H' and not atom2.symbol == 'H':#only bonds between heavy atoms
                        if atom1.sortingLabel < atom2.sortingLabel:
                            bond = molecule.getBond(atom1, atom2)
                            product = CBHSpeciesGenerator().create_cbh1_product(atom1, atom2, bond)
                            spc_list.append(product)
        
        self.map_species_list(self.error_reaction.products, spc_list)
    
    
    def exclude_terminal_atoms(self, molecule, atoms):
        '''
                
        Terminal atoms in the reactant only have atom neighbor so there is
        no overlap of in products. As a result, exclude terminal atoms. 
        '''
        filtered = []
        molecule.sortAtoms()#don't know if this is necessary.
        for atom1 in atoms:
            non_hydrogen_neighbours = []
            for atom2 in atom1.edges:
                    if not atom2.symbol == 'H':
                        non_hydrogen_neighbours.append(atom2)
            if len(non_hydrogen_neighbours) > 1:
                filtered.append(atom1)
                
        return filtered
    
    def account_for_branching(self, molecule, atoms):
        '''
        Atoms with more than two heavy-atom neighbors will 
        not be overcounted once, but more than once.
        
        E.g.
        - a  secondary (non-branched) carbon atom will
        be overcounted once. 
        - a tertiary carbon atom (branched 'once') will
        be overcounted twice.
        
        Generally, the number of reactants per heavy atom 
        that needs to be added to balance the equation
        follows the following formula:
        
        n = no. of heavy atom neighbors - 1
        
        
        '''
        filtered = []
        molecule.sortAtoms()#don't know if this is necessary.
        for atom1 in atoms:
            non_hydrogen_neighbours = []
            for atom2 in atom1.edges:
                    if not atom2.symbol == 'H':
                        non_hydrogen_neighbours.append(atom2)
            filtered.extend([atom1 for i in range(len(non_hydrogen_neighbours)-1)])
                
        return filtered
    def populate_reactants(self):
        '''
        Each pair of adjacent heavy-atom bonds (products in cbh-1)
        shares a common heavy-atom  
        (as reactants in chb-1 as well as products in cbh-0).

        '''
        spc_list = []
        
        #iterate over all heavy-atom bonds:
        molecule = self.spc.molecule[0]
        
        atoms = self.exclude_terminal_atoms(molecule, molecule.atoms)
        atoms = self.exclude_hydrogens(atoms)
        atoms = self.account_for_branching(molecule, atoms)
        
        for atom in atoms:
            reactant = self.createCBH0Product(atom)
            spc_list.append(reactant)
        
        self.map_species_list(self.error_reaction.reactants, spc_list)
    
    
if __name__ == '__main__':
    spc = makeSpeciesFromSMILES('C1CCC(=O)CC1')
    #mol = makeSpecies('C')
    cbh1 = CBH1Reaction(spc=spc)
    cbh1.run()
    rxn =  cbh1.error_reaction
    print rxn.coefficients