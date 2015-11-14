from rmgpy.chemkin import writeKineticsEntry

def write_model(rmg, chemkin_name='reduced_reactions.inp'):
    saveChemkinFile(chemkin_name, rmg.reactionModel.core.species, rmg.reactionModel.core.reactions)


def saveChemkinFile(path, species, reactions, verbose = True):
    

    s ='REACTIONS    KCAL/MOLE   MOLES\n\n'

    for rxn in reactions:
        s +=writeKineticsEntry(rxn, speciesList=species, verbose=verbose)
        s +='\n'
    s +='END\n\n'

    with open(path, 'w') as f:
        f.write(s)