'''
Created on Jul 10, 2014

@author: nickvandewiele
'''
import os
from rmgpy.rmg.main import RMG
from cbh_scheme_generator import *
import logging
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.molecule import Molecule
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.quantity import Quantity
from rmgpy.rmg.input import SMILES, InChI, adjacencyList
from rmgpy import settings
from rmgpy.rmg.main import initializeLog
        
def expand_species_list(species_list, coefficients):
    '''
    Creates a new, expanded list of parse_species, 
    by duplicating the parse_species in the parameter list,
    a number of times equal to the value stored in the
    coefficients dictionary for that particular parse_species.
    
    The duplicated parse_species are references to the original object of 
    the parameter object.
    '''
    expanded_list = []
    for spc in species_list:
        coeff = coefficients[spc.label]
        expanded_list.extend([spc for i in range(coeff)])
    return expanded_list        
################################################################################

class AbstractErrorCancellingReactionParser(object):
    def __init__(self, inputFile):
        '''
        The constructor takes a Species object.
        '''
        self.inputFile = inputFile
        self.rmg = RMG()
        self.error_reaction = None 
        self.speciesDict = {}
        self.global_context = None
        self.local_context = None
    
    def read(self):
        '''
        Reads input on:
        * parse_species whose enthalpy of formation is to be estimated
        * error cancelling reaction: reactants and products of reaction
        * quantum-mechanics method to be used
        * benchmark databases to be used 
        ''' 

        self.parse()
            
        if self.rmg.quantumMechanics:
            self.rmg.quantumMechanics.setDefaultOutputDirectory(self.rmg.outputDirectory)
            self.rmg.reactionModel.quantumMechanics = self.rmg.quantumMechanics    
    
    def read_input_file(self):
        full_path = os.path.abspath(os.path.expandvars(self.inputFile))
        try:
            f = open(full_path)
        except IOError, e:
            logging.error('The input file "{0}" could not be opened.'.format(full_path))
            logging.info('Check that the file exists and that you have read access.')
            raise e
    
        logging.info('Reading input file "{0}"...'.format(full_path))
        try:
            exec f in self.global_context, self.local_context
        except (NameError, TypeError, SyntaxError), e:
            logging.error('The input file "{0}" was invalid:'.format(full_path))
            logging.exception(e)
            raise
        finally:
            f.close()
        logging.info('')  
            
    def initialize_rmg(self):
        self.rmg.reactionModel = CoreEdgeReactionModel()
        self.rmg.initialSpecies = []
        self.rmg.reactionSystems = []
        
    def initialize_context(self):
        self.global_context = { '__builtins__': None }
        self.local_context = {
            '__builtins__': None,
            'True': True,
            'False': False,
            'database': self.parse_database,
            'species': self.parse_species,
            'SMILES': SMILES,
            'InChI': InChI,
            'adjacencyList': adjacencyList,
            'quantumMechanics': self.parse_quantumMechanics,
        }
        
    '''
    TODO
    The parser section of the input file could be much more concise if 
    we could recycle methods from rmgpy.rmg.input.py.
    
    However, that module uses a globally defined rmg object that messes
    up the references to the rmg object here defined. Unless we adapt the methods in
    that module, we cannot use them here.  
    
    '''
    def parse_database(self,
                 thermoLibraries = None,
                 reactionLibraries = None,
                 frequenciesLibraries = None,
                 seedMechanisms = None,
                 kineticsFamilies = 'default',
                 kineticsDepositories = 'default',
                 kineticsEstimator = 'group additivity',
                 ):
        # This function just stores the information about the database to be loaded
        # We don't actually load the database until after we're finished reading
        # the input file
        if isinstance(thermoLibraries, str): thermoLibraries = [thermoLibraries]
        if isinstance(reactionLibraries, str): reactionLibraries = [reactionLibraries]
        if isinstance(seedMechanisms, str): seedMechanisms = [seedMechanisms]
        if isinstance(frequenciesLibraries, str): frequenciesLibraries = [frequenciesLibraries]
        self.rmg.databaseDirectory = settings['database.directory']
        self.rmg.thermoLibraries = thermoLibraries or []
        self.rmg.reactionLibraries = reactionLibraries or []
        self.rmg.seedMechanisms = seedMechanisms or []
        self.rmg.statmechLibraries = frequenciesLibraries or []
        self.rmg.kineticsEstimator = kineticsEstimator
        if kineticsDepositories == 'default':
            self.rmg.kineticsDepositories = ['training']
        elif kineticsDepositories == 'all':
            self.rmg.kineticsDepositories = None
        else:
            assert isinstance(kineticsDepositories,list), "kineticsDepositories should be either 'default', 'all', or a list of names eg. ['training','PrIMe']."
            self.rmg.kineticsDepositories = kineticsDepositories
        if kineticsFamilies in ('default', 'all', 'none'):
            self.rmg.kineticsFamilies = kineticsFamilies
        else:
            assert isinstance(kineticsFamilies,list), "kineticsFamilies should be either 'default', 'all', 'none', or a list of names eg. ['H_Abstraction','R_Recombination'] or ['!Intra_Disproportionation']."
            self.rmg.kineticsFamilies = kineticsFamilies
    
    def parse_species(self, label, structure, reactive=True, unknown=False):
        logging.debug('Found {0} parse_species "{1}" ({2})'.format('reactive' if reactive else 'nonreactive', label, structure.toSMILES()))
        species, isNew = self.rmg.reactionModel.makeNewSpecies(structure, label=label, reactive=reactive)
        props = {
                 'unknown' : unknown
                 }
        species.props = props
        assert isNew, "Species {0} is a duplicate of {1}. Species in input file must be unique".format(label,species.label)
        self.rmg.initialSpecies.append(species)

    def parse_quantumMechanics(self,
                        software,
                        method,
                        fileStore = None,
                        scratchDirectory = None,
                        onlyCyclics = False,
                        maxRadicalNumber = 0,
                        ):
        from rmgpy.qm.main import QMCalculator
        self.rmg.quantumMechanics = QMCalculator()
        self.rmg.quantumMechanics.settings.software = software
        self.rmg.quantumMechanics.settings.method = method
        self.rmg.quantumMechanics.settings.fileStore = fileStore
        self.rmg.quantumMechanics.settings.scratchDirectory = scratchDirectory
        self.rmg.quantumMechanics.settings.onlyCyclics = onlyCyclics
        self.rmg.quantumMechanics.settings.maxRadicalNumber = maxRadicalNumber    

    def run(self):
            if not self.rmg.outputDirectory:
                self.rmg.outputDirectory = os.path.dirname(self.inputFile)
            self.read()

class UserDefinedReactionParser(AbstractErrorCancellingReactionParser):
    def __init__(self, inputFile):
        super(self.__class__, self).__init__(inputFile)
        
    def create_error_cancelling_reaction(self, label, reactants, products, coefficients):
        '''
        Converts the reaction object from the input file into
        a rmg-type reaction by converting the coefficients dict
        into an expanded list of reactants/products.
        '''
        
        for spc in self.rmg.initialSpecies:
            self.speciesDict[spc.label] = spc
            
        #make sure we use references to the same species objects:
        reactants = sorted([self.speciesDict[spec] for spec in reactants])
        products = sorted([self.speciesDict[spec] for spec in products])
        
        expanded_reactants = expand_species_list(reactants, coefficients)
        expanded_products = expand_species_list(products, coefficients)
        
        rxn = Reaction(reactants=expanded_reactants, products=expanded_products)
        self.error_reaction = rxn
        return rxn


    def parse(self):
        """
        Read an error cancelling reaction input file at `path` on disk into the :class:`RMG` object 
        `rmg`.
        """

        self.initialize_rmg()
        self.initialize_context()
        self.local_context['reaction'] = self.create_error_cancelling_reaction
    
        self.read_input_file()
        
    
class UserDefinedSpeciesReactionParser(AbstractErrorCancellingReactionParser):
    def __init__(self, inputFile):
        super(self.__class__, self).__init__(inputFile)
        self.rung = ''
    
    def get_cbh_rxn(self):
        return {
            'cbh0': CBH0Reaction(),
            'cbh1': CBH1Reaction(),
            'cbh2': CBH2Reaction(),
            'cbh3': CBH3Reaction(),
        }[self.rung]
      
    def parse_cbh(self, rung):
        '''
        Converts the reaction object from the input file into
        a rmg-type reaction by converting the coefficients dict
        into an expanded list of reactants/products.
        '''        
        self.rung = rung    
    
    def parse(self):
        """
        Read an error cancelling reaction input file at `path` on disk into the :class:`RMG` object 
        `rmg`.
        """
        self.initialize_rmg()
        self.initialize_context()
        self.local_context['cbh'] = self.parse_cbh
    
        self.read_input_file()
        cbh = self.get_cbh_rxn()
        cbh.spc = self.rmg.initialSpecies[0]
        cbh.run()
        
        #convert ErrorCancellingReaction into a native rmgpy Reaction by expanding list of species
        cbh_rxn = cbh.error_reaction
        
        #put reactants/products in speciesDict:
        for spc in cbh_rxn.reactants:
            self.speciesDict[spc.label] = spc
        for spc in cbh_rxn.products:
            self.speciesDict[spc.label] = spc
            
        #expand initialSpecies with the other reactants and products of the reaction:
        #first re-initialize the list:
        self.rmg.initialSpecies = []
        self.rmg.initialSpecies.extend(self.speciesDict.values())    
        
        reactants = sorted([self.speciesDict[spec.label] for spec in cbh_rxn.reactants])
        products = sorted([self.speciesDict[spec.label] for spec in cbh_rxn.products])
        
        expanded_reactants = expand_species_list(reactants, cbh_rxn.coefficients)
        expanded_products = expand_species_list(products, cbh_rxn.coefficients)
        
        rxn = Reaction(reactants=expanded_reactants, products=expanded_products)
        
        self.error_reaction = rxn

################################################################################


class ErrorCalculator(object):
    def __init__(self, error_reaction, rmg, T):
        self.error_reaction = error_reaction
        self.T = T
        self.rmg = rmg
        self.f_out = ''#output string

    def log_DH298(self, species):
        H = Quantity(species.getEnthalpy(self.T),"J/mol")
        H.units = 'kcal/mol'
        return species.label+' '+ str(H)
        
    def calculate_QM_thermo(self):
        '''
        calculates the thermochemistry of the species
        using QM methods
        '''
        if self.rmg.quantumMechanics:
            self.rmg.quantumMechanics.initialize()
        for species in self.rmg.initialSpecies:
            species.generateThermoData(database=None, quantumMechanics=self.rmg.reactionModel.quantumMechanics)
            logging.info(self.log_DH298(species))
            self.f_out += self.log_DH298(species)+'\n'
        return


    def retrieve_benchmark_thermo(self):
        '''
        Retrieves (hopefully) very accurate values for 
        enthalpy of formation of the reactants and products
        of the error-cancelling reaction.
        '''
        
        self.rmg.loadDatabase()
        
        for spc in self.rmg.initialSpecies:
            if not 'unknown' in spc.props:
                spc.generateThermoData(database=self.rmg.database)
                logging.info(self.log_DH298(spc))
                self.f_out += self.log_DH298(spc)+'\n'
        return
    
    def calculate_improved_enthalpy(self, DHr):
        '''
        Calculates the enthalpy of formation of the specified compound,
        using the reaction enthalpy of the error-cancelling reaction, and
        the benchmark values of the remaining reactants and products 
        of the reaction. 
        '''
        
        value = 0.
        for spc in self.error_reaction.reactants:
            if not 'unknown' in spc.props: 
                value -= spc.getEnthalpy(T)
        
        for spc in self.error_reaction.products:
            value += spc.getEnthalpy(T)
        
        return Quantity(-(DHr-value),"J/mol") 

    def run(self): 
        self.f_out += str(self.error_reaction)+'\n'
        logging.info('QM Enthalpy of Formation: ')
        self.f_out += 'QM Enthalpy of Formation: '+'\n'
        self.calculate_QM_thermo()
        
        DHr = self.error_reaction.getEnthalpyOfReaction(T)
        DHr_kcal = Quantity(DHr, 'J/mol')
        DHr_kcal.units = 'kcal/mol'
        logging.info('Reaction enthalpy: '+str(DHr_kcal))
        self.f_out += 'Reaction enthalpy: '+str(DHr_kcal)+'\n'
        
        logging.info('Benchmark data for enthalpy of formation: ')
        self.f_out += 'Benchmark data for enthalpy of formation: '+'\n'
        self.retrieve_benchmark_thermo()
        
        improved_enthalpy = self.calculate_improved_enthalpy(DHr)
        improved_enthalpy.units = 'kcal/mol'
        unknown_species_label = ''
        for spc in self.error_reaction.reactants:
            if spc.props['unknown']:
                unknown_species_label = spc.label
                break
            
        logging.info('Improved enthalpy of formation of '
                     +unknown_species_label+' '
                     +str(improved_enthalpy)
                     )
        
        self.f_out += 'Improved enthalpy of formation of '+unknown_species_label+': '+str(improved_enthalpy)+'\n'
    
if __name__ == '__main__':
    
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('input', metavar='INPUT', type=str, nargs=1,
        help='Error cancelling reaction input file')
    args = parser.parse_args()
    
    inputFile = os.path.abspath(args.input[0])
    
    level = logging.INFO
    initializeLog(level, 'RMG.log')
    
    #parser = UserDefinedReactionParser(inputFile)
    parser = UserDefinedSpeciesReactionParser(inputFile)
    parser.run()
    
    T = 298
    calc = ErrorCalculator(parser.error_reaction, parser.rmg, T)
    calc.run()
    
    out = 'error_cancelling.out'
    with open(out,'w') as f:
        f.write(calc.f_out)
    
