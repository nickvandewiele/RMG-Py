#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
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

import re
import logging

from rmgpy.molecule import Molecule
from rmgpy.species import Species

# Maximum species label length for Chemkin:
MAX_LENGTH_CHEMKIN_SPECIES_LABEL = 15

# Save species label length for Chemkin:
SAFE_LENGTH_CHEMKIN_SPECIES_LABEL = 10

# Forbidden characters in CHEMKIN species labels:
#The label can only contain alphanumeric characters, hyphens, and underscores
ALLOWED_CHARS = r'[^A-Za-z0-9\-_,\(\)\*]+'

def is_chemkin_compatible(label):
    """
    Checks whether the species label is chemkin compatible by
    verifying a number of conditions:

    - maximum allowed length
    - forbidden characters

    """

    # species length
    if len(label) == 0:
        logging.warning('Species label has length of 0...')
        return False
    elif len(label) > MAX_LENGTH_CHEMKIN_SPECIES_LABEL:
        logging.debug('Species label is longer than 15 characters and will break CHEMKIN 2.0')
        return False
    elif len(label) > SAFE_LENGTH_CHEMKIN_SPECIES_LABEL:
        logging.debug('Species label is longer than 10 characters and may exceed chemkin string limit')

    #forbidden characters
    if not re.search(ALLOWED_CHARS, label):
        return True
    else:
        logging.debug('Species label %s contains a chemkin-incompatible character!', label)
        return False



    return True



def generateSpeciesIdentifier(**kwargs):
    """
    Return a string identifier for the provided augmented inchi that can be used in a
    Chemkin file. Although the Chemkin format allows up to 16 characters for a
    species identifier, this function uses a maximum of 10 to ensure that all
    reaction equations fit in the maximum limit of 52 characters.

    possible kwargs: 
    - species
    - id
    - index
    - user
    - reactive

    """
    
    try:
        species = kwargs['species']
    except KeyError:
        # try if we find an identifier, like an augmented inchi instead
        try:
            identifier = kwargs['id']
            assert 'InChI=1' in identifier, 'Not a valid augmented InChI!'
            species = Species(molecule=[Molecule().fromAugmentedInChI(identifier)])
        except KeyError:
            raise Exception('Could not find an identifier, nor a Species object in {}'.format(kwargs))

    try:
        reactive=kwargs['reactive']
    except KeyError:
        reactive=True

    try:
        index = int(kwargs['index'])
    except KeyError:
        index = -1

    try:
        user_label = kwargs['user']
    except KeyError:
        user_label = ''

    assert len(species.molecule) > 0, 'Species {0} does not contain any molecule objects...'.format(species)

    def get_label(*args, **kwargs): return '{0}'.format(generateLabel(args[0].molecule[0]))
    def get_formula(*args, **kwargs): return '{0}'.format(args[0].molecule[0].getFormula())
    def get_label_index(*args, **kwargs): return '{0}({1:d})'.format(generateLabel(args[0].molecule[0]), args[1])
    def get_formula_index(*args, **kwargs): return '{0}({1:d})'.format(args[0].molecule[0].getFormula(), args[1])
    def get_S_index(*args, **kwargs): return 'S({0:d})'.format(args[1])
    def user_defined(*args, **kwargs): return kwargs['user']

    creators = []
    if not reactive or index == -1:
        creators = [
            get_label, 
            get_formula
            ]
    else:
        creators = [
            get_label_index, 
            get_formula_index, 
            get_S_index
        ]

    if user_label != '':
        creators.insert(0, user_defined)

    for creator in creators:
        proposed_label = creator(species, index, user=user_label)
        if is_chemkin_compatible(proposed_label):
            assert not proposed_label.startswith('('), proposed_label+species.toAdjacencyList()

            return proposed_label

    # If we're here then we just can't come up with a valid Chemkin name
    # for this species, so raise an exception
    raise Exception("Unable to determine valid Chemkin identifier for species {0}.".format(species.toAdjacencyList()))


def generateLabel(mol):
    """
    Use SMILES as default format for label
    However, SMILES can contain slashes (to describe the
    stereochemistry around double bonds); since RMG doesn't 
    distinguish cis and trans isomers, we'll just strip these out
    so that we can use the label in file paths
    """
    smi = mol.toSMILES().replace('/','').replace('\\','')
    label = smi
    
    return label