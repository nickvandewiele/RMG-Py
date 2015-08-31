#!/bin/bash
# Run RMG in parallel with SCOOP on the input files.
python -m scoop ../../rmgpy/reduction/reduction.py minimal/input.py minimal/reduction_input.py minimal/chemkin/chem.inp minimal/chemkin/species_dictionary.txt