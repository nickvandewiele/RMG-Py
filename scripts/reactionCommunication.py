import sys
import os
import numpy as np
import time

import logging
from scoop import logger

from rmgpy import settings
from rmgpy.species import Species
from rmgpy.data.rmg import RMGDatabase
from rmgpy.rmg.react import reactAll

def calculateComm(total, computeTimes, workerCount):
    """

    Calculates the average time spent on communication to/from workers.
    """

    tComm = (workerCount * total - sum(computeTimes)) / workerCount
    logger.info('Average time spent on communication: {:.2f}'.format(tComm))
    logger.info('Contribution of communication to total wall clock time: {:.2f}%'.format(tComm / total * 100))

def setupLogger():

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)

def setupDB():

    db = RMGDatabase()
    db.load(settings['database.directory'], kineticsFamilies='default')

    return db

def loadSpecies(db, library='CHO', useIndex=True):
    """
    Creates a set of molecules from a thermo library.

    Uses the index in the thermo library as the species index.
    Does not 
    """

    spcList = []
    entries = db.thermo.libraries[library].entries

    for label, entry in entries.iteritems():
        index = entry.index if useIndex else -1
        spc = Species(index=index, label=entry.label, molecule = [entry.item])
        spcList.append(spc)
    
    logger.info('No. of species in {} library: {}'.format(library, len(spcList)))
    return spcList

def react(spcList):
    N = len(spcList)
    t1 = time.time()
    rxns, computeTimes, reactionCount = reactAll(spcList, N, np.ones(N), np.ones([N,N]))
    total = time.time() - t1

    logger.info('No. of generated reactions: {}'.format(len(rxns)))
    logger.info('Avg. no. of generated reactions per task: {:.0f}'.format(sum(reactionCount)/len(reactionCount)))
    logger.info('Max. no. of generated reactions per task: {}'.format(max(reactionCount)))

    logger.info('No. of timing results: {}'.format(len(computeTimes)))
    logger.info('Max. time spent in a task: {:.2f}'.format(max(computeTimes)))
    logger.info('Sum of times spent in tasks: {:.2f}'.format(sum(computeTimes)))

    return total, computeTimes    

def main():

    setupLogger()
    db = setupDB()

    spcList = loadSpecies(db)[:100]
    
    fh = logging.FileHandler(filename='timings.log')
    fh.setLevel(logging.INFO)
    logger.addHandler(fh)

    N = len(spcList)
    logger.info(
        'No. of created tasks based on no. of species assuming 1 resonance isomer per species: {}'
        .format(int(N + N*(N+1)/2))
        )
    total, computeTimes = react(spcList)

    logger.info('Total wall clock time: {:.2f}'.format(total))

    workerCount = int(sys.argv[1])
    logger.info('Worker count: {}'.format(workerCount))
    calculateComm(total, computeTimes, workerCount)

    logger.info('Done!')

if __name__ == '__main__':
    main()
