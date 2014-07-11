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
    label='CH4',
    structure=SMILES("C"),
)
species(
    label='C2H6',
    structure=SMILES("CC"),
)

species(
    label='C2H4',
    structure=SMILES("C=C"),
)


reaction(
    label = 'CPD + 5CH4 <=> 3C2H6 + 2C2H4',
    reactants=['CPD','CH4'],
    products=['C2H6', 'C2H4'],
    coefficients={
                  'CPD' :1,
                  'CH4' :5,
                  'C2H6':3,
                  'C2H4':2,
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