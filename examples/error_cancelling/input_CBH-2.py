'''
Isodesmic bond separation error-cancelling reaction for
estimating the enthalpy of formation of 1,3-cyclopentadiene.
'''

database(
    thermoLibraries = ['DFT_QCI_thermo']   
)

species(
    label='CPD',
    structure=SMILES("C1C=CC=C1"),
    unknown='True'
)

species(
    label='C2H6',
    structure=SMILES("CC"),
)

species(
    label='C2H4',
    structure=SMILES("C=C"),
)
species(
    label='C3H8',
    structure=SMILES("CCC"),
)
species(
    label='C3H6',
    structure=SMILES("C=CC"),
)

reaction(
    label = 'CPD + 3C2H6 + 2C2H4 <=> C3H8 + 4C3H4',
    reactants=['CPD','C2H6', 'C2H4'],
    products=['C3H8', 'C3H6'],
    coefficients={
                  'CPD' :1,
                  'C2H6':3,
                  'C2H4':2,
                  'C3H8':1,
                  'C3H6':4
                  }
)

quantumMechanics(
    software='mopac',
    method='pm7',
    fileStore='QMfiles', # relative to where you run it? defaults to inside the output folder.
    scratchDirectory = None, # not currently used
    onlyCyclics = False,
    maxRadicalNumber = 0,
)