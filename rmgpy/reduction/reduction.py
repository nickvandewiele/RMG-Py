#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

#global imports
import copy
import os.path
import numpy as np
import re
import sys
import logging

#global variables
reactions = None


#local imports
from rmgpy.chemkin import getSpeciesIdentifier

from rmgpy.rmg.main import RMG

from rmgpy.reduction.model import ReductionReaction
from rmgpy.reduction.input import load
from rmgpy.reduction.output import write_model
from rmgpy.reduction.rates import isImportant


def simulate_one(reactionModel, atol, rtol, reactionSystem):
    """

    Simulates one reaction system, listener registers results, 
    which are returned at the end.


    The returned data consists of a array of the species names, 
    and the concentration data.

    The concentration data consists of a number of elements for each timestep 
    the solver took to reach the end time of the batch reactor simulation.

    Each element consists of the time and the concentration data of the species at that 
    particular timestep in the order of the species names.

    """

    #register as a listener
    listener = ConcentrationListener()

    coreSpecies = reactionModel.core.species
    regex = r'\([0-9]+\)'#cut of '(one or more digits)'
    species_names = []
    for spc in coreSpecies:
        name = getSpeciesIdentifier(spc)
        name_cutoff = re.split(regex, name)[0]
        species_names.append(name_cutoff)

    listener.species_names = species_names

    reactionSystem.attach(listener)

    pdepNetworks = []
    for source, networks in reactionModel.networkDict.items():
        pdepNetworks.extend(networks)
    
    terminated, obj = reactionSystem.simulate(
        coreSpecies = reactionModel.core.species,
        coreReactions = reactionModel.core.reactions,
        edgeSpecies = reactionModel.edge.species,
        edgeReactions = reactionModel.edge.reactions,
        toleranceKeepInEdge = 0,
        toleranceMoveToCore = 1,
        toleranceInterruptSimulation = 1,
        pdepNetworks = pdepNetworks,
        absoluteTolerance = atol,
        relativeTolerance = rtol,
    ) 

    assert terminated

    #unregister as a listener
    reactionSystem.detach(listener) 

    return listener.species_names, listener.data

def simulate_all(rmg):
    """
    Simulate the RMG job, 
    for each of the simulated reaction systems.

    Each element i of the data corresponds to a reaction system.
    """
    reactionModel = rmg.reactionModel

    data = []

    atol, rtol = rmg.absoluteTolerance, rmg.relativeTolerance
    for _, reactionSystem in enumerate(rmg.reactionSystems):
        data.append(simulate_one(reactionModel, atol, rtol, reactionSystem))

    return data
        

def initialize(wd, rxns):
    global working_dir, reactions
    working_dir = wd
    assert os.path.isdir(working_dir)
    
    #set global variable here such that functions executed in the root worker have access to it.
    
    reactions = [ReductionReaction(rxn) for rxn in rxns]


def find_important_reactions(rmg, tolerance):
    """
    This function:

    - loops over all the species involved in a specific reaction
    - decides whether the specific reaction is important for the species.

    Whenever it is found that a reaction is important for a species, we break
    the species loop, and keep the reaction in the model.


    Returns:
        a list of rxns that can be removed.
    """

    global reactions
    
    # run the simulation, creating concentration profiles for each reaction system defined in input.
    simdata = simulate_all(rmg)

    reduce_reactions = reactions

    """
    Tolerance to decide whether a reaction is unimportant for the formation/destruction of a species

    Tolerance is a floating point value between 0 and 1.

    A high tolerance means that many reactions will be deemed unimportant, and the reduced model will be drastically
    smaller.

    A low tolerance means that few reactions will be deemed unimportant, and the reduced model will only differ from the full
    model by a few reactions.
    """

    CHUNKSIZE = 40
    boolean_array = []
    for chunk in chunks(reduce_reactions,CHUNKSIZE):
        N = len(chunk)
        partial_results = list(
            map(
                assess_reaction, chunk, [rmg.reactionSystems] * N, [tolerance] * N, [simdata] * N
                )
            )
        boolean_array.extend(partial_results)

    important_rxns = []
    for isImport, rxn in zip(boolean_array, reduce_reactions):
        logging.debug('Is rxn {rxn} important? {isImport}'.format(**locals()))
        if isImport:
            important_rxns.append(rxn)


    return [rxn.rmg_reaction for rxn in important_rxns]

def assess_reaction(rxn, reactionSystems, tolerance, data):
    """
    Returns whether the reaction is important or not in the reactions.

    It iterates over the reaction systems, and loads the concentration profile 
    of each reaction system.

    It iterates over a number of samples in profile and 
    evaluates the importance of the reaction at every sample.


    """
    
    logging.debug('Assessing reaction {}'.format(rxn))

    global reactions
    

    # read in the intermediate state variables

    for datum, reactionSystem in zip(data, reactionSystems):    
        T, P = reactionSystem.T.value_si, reactionSystem.P.value_si
        
        species_names, profile = datum

        # take N evenly spaced indices from the table with simulation results:

        """

        Number of time steps between start and end time of the batch reactor simulation at which the importance of 
        reactions should be evaluated.



        The more timesteps, the less chance we have to remove an important reactions, but the more simulations
        need to be carried out.
        """
        
        timesteps = len(profile) / 2
        logging.debug('Evaluating the importance of a reaction at {} time samples.'.format(timesteps))

        assert timesteps <= len(profile)
        indices = map(int, np.linspace(0, len(profile)-1, num = timesteps))
        for index in indices:
            assert profile[index] is not None
            timepoint, coreSpeciesConcentrations = profile[index]

            coreSpeciesConcentrations = {key: float(value) for (key, value) in zip(species_names, coreSpeciesConcentrations)}
            
            for species_i in rxn.reactants:
                if isImportant(rxn, species_i, reactions, 'reactant', tolerance, T, P, coreSpeciesConcentrations):
                    return True

            #only continue if the reaction is not important yet.
            for species_i in rxn.products:
                if isImportant(rxn, species_i, reactions, 'product', tolerance, T, P, coreSpeciesConcentrations):
                    return True

    return False

    
def search_target(target_label, reactionSystem):
    """
    Searches for the Species object in the set of initial species
    that has the same label as the parameter string.
    """

    for k in reactionSystem.initialMoleFractions.keys():
        if k.label == target_label:
            target = k
            break
    assert target is not None, '{} could not be found...'.format(target_label)
    return target


def compute_conversion(target_label, reactionModel, reactionSystem, reactionSystem_index, atol, rtol):
    """
    Computes the conversion of a target molecule by

    - searching the index of the target species in the core species
    of the global reduction variable
    - resetting the reaction system, initialing with empty variables
    - fetching the initial moles variable y0
    - running the simulation at the conditions stored in the reaction system
    - fetching the computed moles variable y
    - computing conversion
    """

    target = search_target(target_label, reactionSystem)
    target_index = reactionModel.core.species.index(target)

    #reset reaction system variables:
    logging.info('No. of rxns in core reactions: {}'.format(len(reactionModel.core.reactions)))
    reactionSystem.initializeModel(\
        reactionModel.core.species, reactionModel.core.reactions,\
        reactionModel.edge.species, reactionModel.edge.reactions, \
        [], atol, rtol)

    #get the initial moles:
    y0 = reactionSystem.y.copy()

    #run the simulation:
    simulate_one(reactionModel, atol, rtol, reactionSystem)

    #compute conversion:
    conv = 1 - (reactionSystem.y[target_index] / y0[target_index])
    return conv

def reduce_compute(tolerance, target_label, reactionModel, rmg, reaction_system_index):
    """
    Reduces the model for the given tolerance and evaluates the 
    target conversion.
    """

    # reduce model with the tolerance specified earlier:
    important_reactions = find_important_reactions(rmg, tolerance)

    original_size = len(reactionModel.core.reactions)

    no_important_reactions = len(important_reactions)
    logging.info('Number of important reactions: {}'.format(no_important_reactions))

    #set the core reactions to the reduced reaction set:
    original_reactions = reactionModel.core.reactions
    rmg.reactionModel.core.reactions = important_reactions

    #re-compute conversion: 
    conversion = compute_conversion(target_label, rmg.reactionModel,\
     rmg.reactionSystems[reaction_system_index], reaction_system_index,\
     rmg.absoluteTolerance, rmg.relativeTolerance)

    #reset the reaction model to its original state:
    rmg.reactionModel.core.reactions = original_reactions

    logging.info('Conversion of reduced model ({} rxns): {:.2f}%'.format(no_important_reactions, conversion * 100))
    return conversion, important_reactions

def optimize_tolerance(target_label, reactionModel, rmg, reaction_system_index, error, orig_conv):
    """
    Increment the trial tolerance from a very low value until the introduced error is greater than the parameter threshold.
    """


    start = 1E-20
    incr = 10
    tolmax = 1

    tol = start
    trial = start

    important_reactions = reactionModel.core.reactions
    
    while True:
        logging.info('Trial tolerance: {trial:.2E}'.format(**locals()))
        Xred, new_important_reactions = reduce_compute(trial, target_label, reactionModel, rmg, reaction_system_index)
        dev = np.abs((Xred - orig_conv) / orig_conv)
        logging.info('Deviation: {dev:.2f}'.format(**locals()))

        if dev > error or trial > tolmax:
            break

        tol = trial
        trial = trial * incr
        important_reactions = new_important_reactions

    if tol == start:
        logging.error('Starting value for tolerance was too high...')

    return tol, important_reactions

class ConcentrationListener(object):
    """Returns the species concentration profiles at each time step."""

    def __init__(self):
        self.species_names = []
        self.data = []

    def update(self, subject):
        self.data.append((subject.t , subject.coreSpeciesConcentrations))


def main():
    
    inputFile, reductionFile, chemkinFile, spc_dict = sys.argv[-4:]

    for f in [inputFile, reductionFile, chemkinFile, spc_dict]:
        assert os.path.isfile(f), 'Could not find {}'.format(f)

    inputDirectory = os.path.abspath(os.path.dirname(inputFile))
    output_directory = inputDirectory

    rmg, target_label, error = load(inputFile, reductionFile, chemkinFile, spc_dict)

    reactionModel = rmg.reactionModel
    initialize(rmg.outputDirectory, reactionModel.core.reactions)

    atol, rtol = rmg.absoluteTolerance, rmg.relativeTolerance
    index = 0
    reactionSystem = rmg.reactionSystems[index]

    print 'Allowed error in target conversion: {0:.0f}%'.format(error * 100)
    
    #compute original target conversion
    Xorig = compute_conversion(target_label, reactionModel, reactionSystem, index,\
     rmg.absoluteTolerance, rmg.relativeTolerance)
    print 'Original target conversion: {0:.0f}%'.format(Xorig * 100)

    # optimize reduction tolerance
    tol, important_reactions = optimize_tolerance(target_label, reactionModel, rmg, index, error, Xorig)
    print 'Optimized tolerance: {:.0E}'.format(tol)

    # plug the important reactions into the RMG object and write:
    rmg.reactionModel.core.reactions = important_reactions
    write_model(rmg)



def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

if __name__ == '__main__':
    main()