'''
Created on Jul 10, 2014

@author: nickvandewiele
'''
import os
from rmgpy.rmg.main import RMG

import logging
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.rmg.input import SMILES, InChI, adjacencyList
from rmgpy import settings
from rmgpy.rmg.main import initializeLog

rmg = None
speciesDict = {}
error_reaction = None
T = 298


################################################################################
'''
TODO
The parser section of the input file could be much more concise if 
we could recycle methods from rmgpy.rmg.input.py.

However, that module uses a globally defined rmg object that messes
up the references to the rmg object here defined. Unless we adapt the methods in
that module, we cannot use them here.  

'''

def parse_database(
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
    rmg.databaseDirectory = settings['database.directory']
    rmg.thermoLibraries = thermoLibraries or []
    rmg.reactionLibraries = reactionLibraries or []
    rmg.seedMechanisms = seedMechanisms or []
    rmg.statmechLibraries = frequenciesLibraries or []
    rmg.kineticsEstimator = kineticsEstimator
    if kineticsDepositories == 'default':
        rmg.kineticsDepositories = ['training']
    elif kineticsDepositories == 'all':
        rmg.kineticsDepositories = None
    else:
        assert isinstance(kineticsDepositories,list), "kineticsDepositories should be either 'default', 'all', or a list of names eg. ['training','PrIMe']."
        rmg.kineticsDepositories = kineticsDepositories
    if kineticsFamilies in ('default', 'all', 'none'):
        rmg.kineticsFamilies = kineticsFamilies
    else:
        assert isinstance(kineticsFamilies,list), "kineticsFamilies should be either 'default', 'all', 'none', or a list of names eg. ['H_Abstraction','R_Recombination'] or ['!Intra_Disproportionation']."
        rmg.kineticsFamilies = kineticsFamilies
    
def parse_species(label, structure, reactive=True, unknown=False):
    logging.debug('Found {0} parse_species "{1}" ({2})'.format('reactive' if reactive else 'nonreactive', label, structure.toSMILES()))
    spec, isNew = rmg.reactionModel.makeNewSpecies(structure, label=label, reactive=reactive)
    spec.props['unknown'] = unknown
    assert isNew, "Species {0} is a duplicate of {1}. Species in input file must be unique".format(label,spec.label)
    rmg.initialSpecies.append(spec)
    speciesDict[label] = spec
        
def parse_quantumMechanics(
                    software,
                    method,
                    fileStore = None,
                    scratchDirectory = None,
                    onlyCyclics = False,
                    maxRadicalNumber = 0,
                    ):
    from rmgpy.qm.main import QMCalculator
    rmg.quantumMechanics = QMCalculator()
    rmg.quantumMechanics.settings.software = software
    rmg.quantumMechanics.settings.method = method
    rmg.quantumMechanics.settings.fileStore = fileStore
    rmg.quantumMechanics.settings.scratchDirectory = scratchDirectory
    rmg.quantumMechanics.settings.onlyCyclics = onlyCyclics
    rmg.quantumMechanics.settings.maxRadicalNumber = maxRadicalNumber    
    
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

def parse_error_cancelling_reaction(label, reactants, products, coefficients):
    '''
    Converts the reaction object from the input file into
    a rmg-type reaction by converting the coefficients dict
    into an expanded list of reactants/products.
    '''
    global speciesDict, error_reaction
    #make sure we use references to the same species objects:
    reactants = sorted([speciesDict[spec] for spec in reactants])
    products = sorted([speciesDict[spec] for spec in products])
    
    expanded_reactants = expand_species_list(reactants, coefficients)
    expanded_products = expand_species_list(products, coefficients)
    
    rxn = Reaction(reactants=expanded_reactants, products=expanded_products)
    error_reaction = rxn
    return rxn
    
def parse_input(path, rmg0):
    """
    Read an error cancelling reaction input file at `path` on disk into the :class:`RMG` object 
    `rmg`.
    """
    global rmg
    full_path = os.path.abspath(os.path.expandvars(path))
    try:
        f = open(full_path)
    except IOError, e:
        logging.error('The input file "{0}" could not be opened.'.format(full_path))
        logging.info('Check that the file exists and that you have read access.')
        raise e

    logging.info('Reading input file "{0}"...'.format(full_path))
    rmg = rmg0
    rmg.reactionModel = CoreEdgeReactionModel()
    rmg.initialSpecies = []
    rmg.reactionSystems = []
    
    global_context = { '__builtins__': None }
    local_context = {
        '__builtins__': None,
        'True': True,
        'False': False,
        'database': parse_database,
        'species': parse_species,
        'SMILES': SMILES,
        'InChI': InChI,
        'adjacencyList': adjacencyList,
        'quantumMechanics': parse_quantumMechanics,
        'reaction': parse_error_cancelling_reaction,
    }

    try:
        exec f in global_context, local_context
    except (NameError, TypeError, SyntaxError), e:
        logging.error('The input file "{0}" was invalid:'.format(full_path))
        logging.exception(e)
        raise
    finally:
        f.close()
    
    logging.info('')    
################################################################################

def read_input(inputFile):
    '''
    Reads input on:
    * parse_species whose enthalpy of formation is to be estimated
    * error cancelling reaction: reactants and products of reaction
    * quantum-mechanics method to be used
    * benchmark databases to be used 
    ''' 
    rmg = RMG()

    if not rmg.outputDirectory:
        rmg.outputDirectory = os.path.dirname(inputFile)
    parse_input(inputFile, rmg)
        
    if rmg.quantumMechanics:
        rmg.quantumMechanics.setDefaultOutputDirectory(rmg.outputDirectory)
        rmg.reactionModel.quantumMechanics = rmg.quantumMechanics
            
    return rmg

def log_DH298(species):
    logging.info(species.label+' '+ str(species.getEnthalpy(T)))
    
def calculate_QM_thermo():
    '''
    calculates the thermochemistry of the species
    using QM methods
    '''
    global rmg,speciesDict
    if rmg.quantumMechanics:
        rmg.quantumMechanics.initialize()
    for species in rmg.initialSpecies:
        species.generateThermoData(database=None, quantumMechanics=rmg.reactionModel.quantumMechanics)
        log_DH298(species)
    return


def retrieve_benchmark_thermo():
    '''
    Retrieves (hopefully) very accurate values for 
    enthalpy of formation of the reactants and products
    of the error-cancelling reaction.
    '''
    
    rmg.loadDatabase()
    
    for spc in rmg.initialSpecies:
        if not spc.props['unknown']:
            spc.generateThermoData(database=rmg.database)
            log_DH298(spc)
    return

def calculate_improved_enthalpy(DHr):
    '''
    Calculates the enthalpy of formation of the specified compound,
    using the reaction enthalpy of the error-cancelling reaction, and
    the benchmark values of the remaining reactants and products 
    of the reaction. 
    '''
    global error_reaction
    
    value = 0.
    for spc in error_reaction.reactants:
        if not spc.props['unknown']: 
            value -= spc.getEnthalpy(T)
    
    for spc in error_reaction.products:
        value += spc.getEnthalpy(T)
    
    return -(DHr-value)

def run(inputFile):
    global error_reaction
    
    read_input(inputFile)
    logging.info('QM Enthalpy of Formation: ')
    calculate_QM_thermo()
    
    DHr = error_reaction.getEnthalpyOfReaction(T)
    logging.info('Reaction enthalpy: '+str(DHr))
    
    logging.info('Benchmark data for enthalpy of formation: ')
    retrieve_benchmark_thermo()
    
    improved_enthalpy = calculate_improved_enthalpy(DHr)
    
    unknown_species_label = ''
    for spc in error_reaction.reactants:
        if spc.props['unknown']:
            unknown_species_label = spc.label
            break
        
    logging.info('Improved enthalpy of formation of '
                 +unknown_species_label+' '
                 +str(improved_enthalpy)
                 )
    
if __name__ == '__main__':
    
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('input', metavar='INPUT', type=str, nargs=1,
        help='Error cancelling reaction input file')
    args = parser.parse_args()
    
    inputFile = os.path.abspath(args.input[0])
    
    level = logging.INFO
    initializeLog(level, 'RMG.log')
    
    run(inputFile)