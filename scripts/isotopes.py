#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script accepts two input files, one input file (e.g. input.py) with the RMG-Py model to generate, NOT
containing any non-normal isotopomers, and one input file (e.g. input_isotope.py) for the model to be 
generated starting from the isotopomers generated from the core species of the first model.
"""

import argparse
import logging
import numpy as np
import os
import os.path

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
# General syntax to import a library but no functions: 
##import (library) as (give the library a nickname/alias)
import matplotlib.pyplot as plt

from rmgpy.tools.isotopes import run

################################################################################


def parseCommandLineArguments():
    """
    Parse the command-line arguments being passed to RMG-Py. This uses the
    :mod:`argparse` module, which ensures that the command-line arguments are
    sensible, parses them, and returns them.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='RMG input file')
    parser.add_argument('isotopeinput', help='RMG input file for the isotope model')
    parser.add_argument('output', help='Output folder')
    args = parser.parse_args()
    
    return args

def findSMILES(label, spcList):
    """
    Searches through the list of species and returns the SMILES of the species
    whose label + index is equal to the parameter string.
    """

    match = None
    for spc in spcList:
        if spc.index == -1:
            if '{}'.format(spc.label) == label:
                match = spc
                break
        else:            
            if '{}({})'.format(spc.label, spc.index) == label:
                match = spc
                break

    mol = match.molecule[0]
    if mol.getNumAtoms() - mol.getNumAtoms('H') == 0: return None

    if match.reactive:
        print 'label: {}'.format(label)
        return mol.toSMILES() 
    else: return None

    raise Exception('Could not find the label {} in the list of species...'.format(label))

def plot(probs, spcList, timeData, outputDirectory):
    """
    
    """

    for df in probs:
        spcLabel = df.columns[0]
        smi = findSMILES(spcLabel, spcList)
        if smi:
            plt.figure()
            fig, ax = plt.subplots()
            for label in df:
                plt.plot(timeData, df[label], lw=1)
            
            plt.title(smi)
            plt.legend(loc='best')
            plt.ylim([0,1])
            plt.xlabel('Time (s)')
            plt.ylabel('Probability (-)')
            plt.grid(b=True, which='major')
            plt.yticks(np.arange(0, 1.1, 0.1))


            print 'Saving figure: {}'.format(os.path.abspath(smi)+'.png')
            plt.savefig(os.path.join(outputDirectory, smi+'.png'))

def main():

    args = parseCommandLineArguments()

    inputFile = args.input
    inputIsoFile = args.isotopeinput
    outputdir = os.path.abspath(args.output)

    probs, spcs, spcdata = run(inputFile, inputIsoFile, outputdir)
    
    figs = os.path.join(outputdir, 'figs')
    os.mkdir(figs)
    print ('Creating figures in {}.'.format(figs))
    plot(probs, spcs, spcdata['Time (s)'], figs)

if __name__ == '__main__':
    main()