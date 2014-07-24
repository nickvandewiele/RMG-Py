'''
Created on Jul 21, 2014

@author: nickvandewiele
'''

from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import CoreEdgeReactionModel#TODO ideally, we don't want to use CoreEdgeReactionModel

#####################################################################
#   Utiliy Methods
#####################################################################
def exclude_hydrogens(atoms):
        return [atom for atom in atoms if not atom.symbol == 'H']

def is_terminal_atom(atom):
    '''
    This atom has less than 2 heavy atoms attached
    '''
    neighbors = exclude_hydrogens([atom1 for atom1 in atom.edges])
    return len(neighbors) < 2

def is_terminal_bond(atom1, atom2):
    return is_terminal_atom(atom1) or is_terminal_atom(atom2)

def is_connected_to_terminal_bond(atom):
        '''
        If any of its neighboring atoms is a terminal atom,
        this method returns True.
    
        '''
        atoms = exclude_hydrogens(atom.edges)
        for atom2 in atoms:
            if is_terminal_atom(atom2):
                return True
            
        return False

def make_species_from_adjacencyList(adjList):
    mol = Molecule().fromAdjacencyList(adjList, saturateH=True)#sature with hydrogens
    product = make_species_from_molecule(mol)
    return product

def make_species_from_molecule(mol, label=''):
    spc, isNew = CoreEdgeReactionModel().makeNewSpecies(mol, label=label)
    return spc

def make_species_from_SMILES(smi):
    mol = Molecule().fromSMILES(smi)
    return make_species_from_molecule(mol, label=smi)

def get_all_bonds(molecule):
    bonds = []
    for atom1 in molecule.vertices:
            for atom2 in atom1.edges:
                if atom1.sortingLabel < atom2.sortingLabel:
                    bonds.append([atom1, atom2])
    return bonds

def exclude_hydrogen_bonds(bonds):
    filtered = []
    for bond in bonds:
        contains_h = False
        for atom in bond:
            if atom.symbol == 'H':
                contains_h = True
                break
        if not contains_h:
            filtered.append(bond)
            
    return filtered

def exclude_terminal_bonds(bonds):
    filtered = []
    for b in bonds:
        atom1, atom2 = b[0], b[1]
        if is_terminal_bond(atom1, atom2):# do not include terminal atoms
            pass
        else:
            filtered.append(b)
    return filtered            


#####################################################################


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
    '''
    TODO Ideally we want to create species WITHOUT having the create
    a de-tour by writing an adjacency list that is then parsed again.
    
    However, because of the difficulty of creating molecules in RMG-Py
    i.e., there are no wrapper methods like 'molecule.addBond(atom1, atom2) that
    will do the 'dirty' accounting (e.g. index updating, atom typing, etc...)
    for me, I resort to creating adjacency lists instead.
     
    '''
    def __init__(self):
        pass
    
    def write_electron_info(self, atom):
        '''
        writes the part:
        C 1 0 
        
        in:
        
        1 C 1 0 {2,S}{...}
        '''
        return ' '.join([atom.symbol,str(atom.radicalElectrons),str(atom.lonePairs)])
                        
    def determine_index(self, atom, index_map, running_index):
        '''
        This method looks up in the parameter map 'index_map' whether
        the parameter atom 'atom' has already an index assigned to it.
        
        If so, the return that value, and don't change the value of the
        parameter 'running_index'.
        
        If that atom was not yet encountered, then we will assign it 
        a new value, equal to the value of 'running_index'+1.
        
        Eventually, we would like to return the retrieved (or generated)
        index of the parameter atom AND
        the updated value of the 'running_index'.
        '''
        if atom in index_map.keys():
            index =  index_map[atom]
        else:
            running_index += 1#update
            index = running_index#assign the index to the updated value
            index_map[atom] = index#put that newly encountered atom in the map
            
        return index, running_index
    
    def write_connectivity_cbh3_product(self, central_atom, adjacent_to_central_atom, indices, running_index, molecule, lines):
        '''
        This method takes a central atom and will do two things to build a part of the adjacency list:
        
        1) it will create a line with the connectivity info of the central atom
        E.g.:
        
        1 C 0 0 {2,S} {3,S}
        
        2) it will create a line with the connectivity info the neighbors of that central atom,
        excluding that one atom that is 'special', i.e. the variable 'adjacent_to_central_atom'.
        E.g.:
        
        3 C 0 0 {1,S}
        
        This method has to retrieve the indices of the neighboring atoms of the central atoms,
        and delegates the retrieval of the index of that atom to another method.
        
        
        '''
        #find the neighbors of the central atom, except that one 'special' neighbor:
        other_neighbors = [atom for atom in central_atom.edges if not atom == adjacent_to_central_atom]
        other_neighbors = exclude_hydrogens(other_neighbors)

        connectivity = '{'+str(indices[adjacent_to_central_atom])+','+molecule.getBond(central_atom, adjacent_to_central_atom).order+'}'
        for neigh in other_neighbors:
            order = molecule.getBond(central_atom, neigh).order
            index, running_index = self.determine_index(neigh, indices, running_index)
            connectivity = ' '.join([connectivity, '{'+str(index)+','+order+'}'])
            e_info = self.write_electron_info(neigh)
            line = ' '.join([str(index),e_info,'{'+str(indices[central_atom])+','+order+'}'])
            lines.append(line)
        
        #create the line of the central atom:
        e_info = self.write_electron_info(central_atom)
        line_central_atom = ' '.join([str(indices[central_atom]),e_info, connectivity])
        
        #inserts line of central atom at the top of the adjacency list, or just the 2nd line:
        lines.insert(indices[central_atom]-1, line_central_atom) 
    
        return running_index
    
    def createAdjacencyList_cbh1_product(self, atom1, atom2, bond):
        order = bond.order
        e_info = self.write_electron_info(atom1)
        first_line = '1 '+e_info+' {2,'+order+'}'
        e_info = self.write_electron_info(atom2)
        second_line = '2 '+e_info+' {1,'+order+'}'
        adjList = first_line + '\n' + second_line
        
        return adjList
    
    def createAdjacencyList_cbh2_product(self, atom, neighbors, molecule):
        lines = []
        e_info = self.write_electron_info(atom)
        first_line = ' '.join(['1',e_info])
        connectivity = ''
        for i, neigh in enumerate(neighbors):
            order = molecule.getBond(atom, neigh).order
            index_neigh = str(i+2)#neighbours start at index 2
            connectivity += '{'+index_neigh+','+order+'} '
            e_info = self.write_electron_info(neigh)
            line = ' '.join([index_neigh,e_info,'{1,'+order+'}'])
            lines.append(line)
            
        first_line += ' '+connectivity
        lines.insert(0, first_line)
        
        return '\n'.join(lines)

    
    def createAdjacencyList_cbh3_product(self, atom1, atom2, molecule):
        lines = []
    
        
        '''
        the variable 'indices' is a map that keeps track of the assigned indices to the 
        2 central atoms, and their neighbors.
        
        E.g.:
        
        the neighbor of atom1 will be assigned index '3', which should be stored,
        in case atom2 is also connected to that same atom.
        '''
        #we initiate the map with the values for the indices of the atoms in the 
        #new adjacency list we already know:
        indices = {
                     atom1 : 1,
                     atom2 : 2
                     }
        
        #initialize the running index:
        running_index = 2 #the two central atoms will have index '1' and '2'
        
        running_index = self.write_connectivity_cbh3_product(atom1, atom2, indices, running_index, molecule, lines)
        running_index = self.write_connectivity_cbh3_product(atom2, atom1, indices, running_index, molecule, lines)
        
        return '\n'.join(lines)
    
    def create_cbh0_product(self, atom):
        return make_species_from_SMILES(atom.symbol)
        
    def create_cbh1_product(self, atom1, atom2, bond):
        adjList = self.createAdjacencyList_cbh1_product(atom1, atom2, bond)#possibly a de-tour by creating adjList
        product = make_species_from_adjacencyList(adjList)
        return product 
    
    def create_cbh2_product(self, atom, neighbors, molecule):
        adjList = self.createAdjacencyList_cbh2_product(atom, neighbors, molecule)#possibly a de-tour by creating adjList
        product = make_species_from_adjacencyList(adjList)
        return product
    
    def create_cbh3_product(self, atom1, atom2, molecule):
        adjList = self.createAdjacencyList_cbh3_product(atom1, atom2, molecule)#possibly a de-tour by creating adjList
        product = make_species_from_adjacencyList(adjList)
        return product
    
class Abstract_CBH_Reaction(object):
    '''
    Superclass to all implementations of the 
    Connectivity-based Hierarchy error-canceling reaction classes.
    
    '''
    def __init__(self, spc=None):
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


class CBH0Reaction(Abstract_CBH_Reaction):
    '''
    Creates rung '0' of the CBH method for the creation
    of an potential error-canceling reaction.
    '''  
    def __init__(self, spc=None):
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
        heavy_atoms = exclude_hydrogens(molecule.atoms)
        for atom in heavy_atoms:#iterate over all atoms!
            product = make_species_from_SMILES(atom.symbol)
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
        
        bonds = get_all_bonds(molecule)
        
        bonds = exclude_hydrogen_bonds(bonds)
        
        for b in bonds:
            bond = molecule.getBond(b[0], b[1])
            balancing_dihydrogen += bond_orders[bond.order]
    
        h2 = make_species_from_SMILES('[H][H]')
        self.error_reaction.reactants.append(h2)
        self.error_reaction.coefficients[h2.label] = balancing_dihydrogen
    

class CBH1Reaction(Abstract_CBH_Reaction):
    '''
    Creates rung '1' of the CBH method for the creation
    of an potential error-canceling reaction.
    
    Corresponds to the Pople's isodesmic bond separation scheme.
    '''  
    def __init__(self, spc=None):
        super(self.__class__, self).__init__(spc)
        
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
            neighbours = exclude_hydrogens([atom2 for atom2 in atom1.edges])
            filtered.extend([atom1 for _ in range(len(neighbours)-1)])
                
        return filtered
        
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
    
    def populate_reactants(self):
        '''
        Each pair of adjacent heavy-atom bonds (products in cbh-1)
        shares a common heavy-atom  
        (as reactants in chb-1 as well as products in cbh-0).

        '''
        spc_list = []
        
        #iterate over all heavy-atom bonds:
        molecule = self.spc.molecule[0]
        
        atoms = exclude_hydrogens(molecule.atoms)
        atoms = self.account_for_branching(molecule, atoms)
        
        for atom in atoms:
            reactant = CBHSpeciesGenerator().create_cbh0_product(atom)
            spc_list.append(reactant)
        
        self.map_species_list(self.error_reaction.reactants, spc_list)

class CBH2Reaction(Abstract_CBH_Reaction):
    '''
    Creates rung '2' of the CBH method for the creation
    of an potential error-canceling reaction.
    
    Corresponds to the simplest hypohomodesmotic reaction scheme
    developed by Wheeler et al. 
    
    Another, more appropriate, illuminative name is isoatomic scheme.
    '''  
    def __init__(self, spc=None):
        super(self.__class__, self).__init__(spc)
    
    def populate_products(self):  
        '''
        Atom centered method. 
        Preserve atom environment of atoms.
        
        Exclude terminal atoms or hydrogens.
        '''

        spc_list = []
        
        molecule = self.spc.molecule[0]
        atoms = exclude_hydrogens(molecule.atoms)
        
        #iterate over all heavy-atom atoms:
        for atom in atoms:
            neighbors = [atom2 for atom2 in atom.edges]
            neighbors = exclude_hydrogens(neighbors)
            if not is_terminal_atom(atom):#don't include terminal atoms
                product = CBHSpeciesGenerator().create_cbh2_product(atom, neighbors, molecule)
                spc_list.append(product)
        
        self.map_species_list(self.error_reaction.products, spc_list)

        
    def populate_reactants(self):
        '''
        Each pair of adjacent heavy-atom atoms
        preserving their immediate environment (products in cbh-2) shares 
        a common heavy-atom bond 
        (reactants in cbh2 as well as products in cbh-1)
        
        Exclude terminal bonds, or bonds with hydrogen.

        '''
        spc_list = []
        
        molecule = self.spc.molecule[0]
        
        #iterate over all heavy-atom bonds:
        molecule.sortAtoms()
        
        bonds = get_all_bonds(molecule)
        
        bonds = exclude_hydrogen_bonds(bonds)
        
        bonds = exclude_terminal_bonds(bonds)
        for b in bonds:
            atom1, atom2 = b[0], b[1]
            bond = molecule.getBond(atom1, atom2)
            reactant = CBHSpeciesGenerator().create_cbh1_product(atom1, atom2, bond)
            spc_list.append(reactant)
        
        self.map_species_list(self.error_reaction.reactants, spc_list)

class CBH3Reaction(Abstract_CBH_Reaction):
    '''
    Creates rung '3' of the CBH method for the creation
    of an potential error-canceling reaction.
    
    Corresponds to the simplest hyperhomodesmotic reaction scheme
    developed by Wheeler et al., and preserve the immediate connectivity of all bonds in the molecule,
    i.e. every heavy-atom bond is extracted maintaining its immediate connectivity.
    '''  
    def __init__(self, spc=None):
        super(self.__class__, self).__init__(spc)
    
    def account_for_branching(self, molecule, atoms):
        '''
        
        Generally, the number of reactants per heavy atom 
        that needs to be added to balance the equation
        follows the following formula:
        
        n = no. of heavy atom neighbors
            - 1 (if atom belongs to terminal bond)
            - 1
        
        
        '''
        filtered = []
        molecule.sortAtoms()#don't know if this is necessary.
        for atom1 in atoms:
            neighbors = exclude_hydrogens([atom2 for atom2 in atom1.edges])
            
            n = len(neighbors)-1#branching
            if is_connected_to_terminal_bond(atom1):
                n -= 1
                
            filtered.extend([atom1 for _ in range(n)])
                
        return filtered
    
    
    def populate_products(self): 
        '''
        
        Bond centered method.
        Preserve bond environment of each bond.
        
        Exclude terminal bonds, or bonds with hydrogen.
        ''' 
        spc_list = []

        #iterate over all heavy-atom bonds:
        molecule = self.spc.molecule[0]
        
        #iterate over all unique non-hydrogen bonds of the Molecule:
        molecule.sortAtoms()#don't know if this is necessary.
        bonds = get_all_bonds(molecule)
        
        bonds = exclude_hydrogen_bonds(bonds)
        
        bonds = exclude_terminal_bonds(bonds)
        
        for b in bonds:
            atom1, atom2 = b[0], b[1]
            product = CBHSpeciesGenerator().create_cbh3_product(atom1, atom2, molecule)
            spc_list.append(product)
        
        self.map_species_list(self.error_reaction.products, spc_list)
        
    def populate_reactants(self):
        '''
        
        Exclude hydrogens and account for branching.
        '''
        spc_list = []
        
        #iterate over all heavy-atom bonds:
        molecule = self.spc.molecule[0]
        
        atoms = exclude_hydrogens(molecule.atoms)
        atoms = self.account_for_branching(molecule, atoms)
        
        for atom in atoms:
            neighbors = exclude_hydrogens([atom2 for atom2 in atom.edges])
            reactant = CBHSpeciesGenerator().create_cbh2_product(atom, neighbors, molecule)
            spc_list.append(reactant)
        
        self.map_species_list(self.error_reaction.reactants, spc_list)
         
        
if __name__ == '__main__':
    spc = make_species_from_SMILES('C1C=CC=C1')
    #mol = makeSpecies('C')
    cbh = CBH0Reaction(spc=spc)
    cbh.run()
    rxn =  cbh.error_reaction
    print rxn.coefficients