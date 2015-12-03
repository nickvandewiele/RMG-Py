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

import sys
import os.path

from input import load
from output import write_model
from optimization import optimize
from reduction import compute_observables, initialize

from rmgpy.scoop_framework.util import logger

def main():
    level = 20#10 : debug, 20: info
    initializeLog(level)
    inputFile, reductionFile, chemkinFile, spc_dict = sys.argv[-4:]

    for f in [inputFile, reductionFile, chemkinFile, spc_dict]:
        assert os.path.isfile(f), 'Could not find {}'.format(f)

    inputDirectory = os.path.abspath(os.path.dirname(inputFile))
    output_directory = inputDirectory

    rmg, targets, error = load(inputFile, reductionFile, chemkinFile, spc_dict)
    logger.info('Allowed error in target observables: {0:.0f}%'.format(error * 100))

    reactionModel = rmg.reactionModel
    initialize(rmg.outputDirectory, reactionModel.core.reactions)

    atol, rtol = rmg.absoluteTolerance, rmg.relativeTolerance
    index = 0
    reactionSystem = rmg.reactionSystems[index]    
    
    #compute original target observables
    observables = compute_observables(targets, reactionModel, reactionSystem, \
     rmg.absoluteTolerance, rmg.relativeTolerance)

    logger.info('Observables of original model:'.format())
    for target, observable in zip(targets, observables):
        logger.info('{}: {:.2f}%'.format(target, observable * 100))

    # optimize reduction tolerance
    tol, important_reactions = optimize(targets, reactionModel, rmg, index, error, observables)
    logger.info('Optimized tolerance: {:.0E}'.format(10**tol))
    logger.info('Number of reactions in optimized reduced model : {}'.format(len(important_reactions)))

    # plug the important reactions into the RMG object and write:
    rmg.reactionModel.core.reactions = important_reactions
    write_model(rmg)

def initializeLog(level):
    """
    Set up a logger for reduction to use to print output to stdout. The
    `level` parameter is an integer specifying the amount of log text seen
    at the console; the levels correspond to those of the :data:`logging` module.
    """
    # Create logger
    logger.setLevel(level)


if __name__ == '__main__':
    main()