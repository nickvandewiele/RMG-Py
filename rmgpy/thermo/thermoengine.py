
import numpy
import math

import logging as logging
    

from rmgpy.thermo.async import AbstractEngine
import rmgpy.constants as constants
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.statmech import Conformer
from rmgpy.thermo import Wilhoit, NASA, ThermoData
import rmgpy.data.rmg

thermo_engine = None

def init():
    """
    Creates a global variable thermo_engine 
    and initializes it.
    """
    global thermo_engine
    thermo_engine = ThermoEngine()

def isValid(identifier):
    """Check if the parameter string is a valid identifier"""
    return 'InChI=' in identifier

def parse_inchi(aug_inchi):
    """Convert the identifier to a dictionary key."""
    return aug_inchi.split('/', 1)[1]

def processThermoData(spc, thermo0, thermoClass=NASA):
    """
    Converts via Wilhoit into required `thermoClass` and sets `E0`.
    
    Resulting thermo is stored (`spc.thermo`) and returned.
    """

    database = rmgpy.data.rmg.database
    solvationdatabase = database.solvation

    # Always convert to Wilhoit so we can compute E0
    if isinstance(thermo0, Wilhoit):
        wilhoit = thermo0
    elif isinstance(thermo0, ThermoData):
        Tdata = thermo0._Tdata.value_si
        Cpdata = thermo0._Cpdata.value_si
        H298 = thermo0._H298.value_si
        S298 = thermo0._S298.value_si
        Cp0 = thermo0._Cp0.value_si
        CpInf = thermo0._CpInf.value_si
        wilhoit = Wilhoit().fitToDataForConstantB(Tdata, Cpdata, Cp0, CpInf, H298, S298, B=1000.0)
    else:
        Cp0 = spc.calculateCp0()
        CpInf = spc.calculateCpInf()
        wilhoit = thermo0.toWilhoit(Cp0=Cp0, CpInf=CpInf)
    wilhoit.comment = thermo0.comment

    # Add on solvation correction
    if LiquidModel.solventData and not "Liquid thermo library" in thermo0.comment:
        #logging.info("Making solvent correction for {0}".format(LiquidModel.solventName))
        soluteData = solvationdatabase.getSoluteData(spc)
        solvation_correction = solvationdatabase.getSolvationCorrection(soluteData, LiquidModel.solventData)
        # correction is added to the entropy and enthalpy
        wilhoit.S0.value_si = (wilhoit.S0.value_si + solvation_correction.entropy)
        wilhoit.H0.value_si = (wilhoit.H0.value_si + solvation_correction.enthalpy)
        
    # Compute E0 by extrapolation to 0 K
    if spc.conformer is None:
        spc.conformer = Conformer()
    spc.conformer.E0 = wilhoit.E0
    
    # Convert to desired thermo class
    if thermoClass is Wilhoit:
        spc.thermo = wilhoit
    elif thermoClass is NASA:
        if LiquidModel.solventData:
            #if liquid phase simulation keep the nasa polynomial if it comes from a liquid phase thermoLibrary. Otherwise convert wilhoit to NASA
            if "Liquid thermo library" in thermo0.comment and isinstance(thermo0, NASA):
                spc.thermo = thermo0
                if spc.thermo.E0 is None:
                    spc.thermo.E0 = wilhoit.E0
            else:
                spc.thermo = wilhoit.toNASA(Tmin=100.0, Tmax=5000.0, Tint=1000.0)
        else: 
            #gas phase with species matching thermo library keep the NASA from library or convert if group additivity
            if "Thermo library" in thermo0.comment and isinstance(thermo0,NASA):
                spc.thermo=thermo0
                if spc.thermo.E0 is None:
                    spc.thermo.E0 = wilhoit.E0
            else:
                spc.thermo = wilhoit.toNASA(Tmin=100.0, Tmax=5000.0, Tint=1000.0)
    else:
        raise Exception('thermoClass neither NASA nor Wilhoit.  Cannot process thermo data.')
    
    if spc.thermo.__class__ != thermo0.__class__:
        # Compute RMS error of overall transformation
        Tlist = numpy.array([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0], numpy.float64)
        err = 0.0
        for T in Tlist:
            err += (spc.thermo.getHeatCapacity(T) - thermo0.getHeatCapacity(T))**2
        err = math.sqrt(err/len(Tlist))/constants.R
        # logging.log(logging.WARNING if err > 0.2 else 0, 'Average RMS error in heat capacity fit to {0} = {1:g}*R'.format(spc, err))

    return spc.thermo
    

def generateThermoData(spc, thermoClass=NASA):
    """
    Generates thermo data, first checking Libraries, then using either QM or Database.
    
    If quantumMechanics is not None, it is asked to calculate the thermo.
    Failing that, the database is used.
    
    The database generates the thermo data for each structure (resonance isomer),
    picks that with lowest H298 value.
    
    It then calls :meth:`processThermoData`, to convert (via Wilhoit) to NASA
    and set the E0.
    
    Result stored in `spc.thermo` and returned.
    """
    

    database = rmgpy.data.rmg.database
    quantumMechanics = None #rmg.quantumMechanics TODO this should be reset for production purposes

    thermo0 = None
    
    thermodatabase = database.thermo
    assert thermodatabase is not None

    thermo0 = thermodatabase.getThermoDataFromLibraries(spc)
    
    if thermo0 is not None:
        logging.info("Found thermo for {0} in {1}".format(spc.label,thermo0[0].comment.lower()))
        assert len(thermo0) == 3, "thermo0 should be a tuple at this point: (thermoData, library, entry)"
        thermo0 = thermo0[0]
        
    elif quantumMechanics:
        original_molecule = spc.molecule[0]
        if quantumMechanics.settings.onlyCyclics and not original_molecule.isCyclic():
            pass
        else: # try a QM calculation
            if original_molecule.getRadicalCount() > quantumMechanics.settings.maxRadicalNumber:
                # Too many radicals for direct calculation: use HBI.
                logging.info("{0} radicals on {1} exceeds limit of {2}. Using HBI method.".format(
                    original_molecule.getRadicalCount(),
                    spc.label,
                    quantumMechanics.settings.maxRadicalNumber,
                    ))
                
                # Need to estimate thermo via each resonance isomer
                thermo = []
                for molecule in spc.molecule:
                    molecule.clearLabeledAtoms()
                    # Try to see if the saturated molecule can be found in the libraries
                    tdata = thermodatabase.estimateRadicalThermoViaHBI(molecule, thermodatabase.getThermoDataFromLibraries)
                    priority = 1
                    if tdata is None:
                        # Then attempt quantum mechanics job on the saturated molecule
                        tdata = thermodatabase.estimateRadicalThermoViaHBI(molecule, quantumMechanics.getThermoData)
                        priority = 2
                    if tdata is None:
                        # Fall back to group additivity
                        tdata = thermodatabase.estimateThermoViaGroupAdditivity(molecule)
                        priority = 3
                    
                    thermo.append((priority, tdata.getEnthalpy(298.), molecule, tdata))
                
                if len(thermo) > 1:
                    # Sort thermo first by the priority, then by the most stable H298 value
                    thermo = sorted(thermo, key=lambda x: (x[0], x[1])) 
                    # Save resonance isomers reordered by their thermo
                    spc.molecule = [item[2] for item in thermo]
                    original_molecule = spc.molecule[0]
                thermo0 = thermo[0][3] 
                
                # If priority == 2
                if thermo[0][0] == 2:
                    # Write the QM molecule thermo to a library so that can be used in future RMG jobs.  (Do this only if it came from a QM calculation)
                    quantumMechanics.database.loadEntry(index = len(quantumMechanics.database.entries) + 1,
                                                    label = original_molecule.toSMILES() + '_({0})'.format(_multiplicity_labels[original_molecule.multiplicity]),
                                                    molecule = original_molecule.toAdjacencyList(),
                                                    thermo = thermo0,
                                                    shortDesc = thermo0.comment
                                                    
                                                    )                    
#                    # For writing thermodata HBI check for QM molecules
#                    with open('thermoHBIcheck.txt','a') as f:
#                        f.write('// {0!r}\n'.format(thermo0).replace('),','),\n//           '))
#                        f.write('{0}\n'.format(original_molecule.toSMILES()))
#                        f.write('{0}\n\n'.format(original_molecule.toAdjacencyList(removeH=False)))

            else: # Not too many radicals: do a direct calculation.
                thermo0 = quantumMechanics.getThermoData(original_molecule) # returns None if it fails
            
                if thermo0 is not None:
                    # Write the QM molecule thermo to a library so that can be used in future RMG jobs.
                    quantumMechanics.database.loadEntry(index = len(quantumMechanics.database.entries) + 1,
                                                    label = original_molecule.toSMILES() + '_({0})'.format(_multiplicity_labels[original_molecule.multiplicity]),
                                                    molecule = original_molecule.toAdjacencyList(),
                                                    thermo = thermo0,
                                                    shortDesc = thermo0.comment
                                                    )                    
    if thermo0 is None:
        # Use group additivity methods to determine thermo for molecule (or if QM fails completely)
        original_molecule = spc.molecule[0]
        if original_molecule.getRadicalCount() > 0:
            # Molecule is a radical, use the HBI method
            thermo = []
            for molecule in spc.molecule:
                molecule.clearLabeledAtoms()
                # First see if the saturated molecule is in the libaries
                tdata = thermodatabase.estimateRadicalThermoViaHBI(molecule, thermodatabase.getThermoDataFromLibraries)
                priority = 1
                if tdata is None:
                    # Otherwise use normal group additivity to obtain the thermo for the molecule
                    tdata = thermodatabase.estimateThermoViaGroupAdditivity(molecule)
                    priority = 2
                thermo.append((priority, tdata.getEnthalpy(298.), molecule, tdata))
            
            if len(thermo) > 1:
                # Sort thermo first by the priority, then by the most stable H298 value
                thermo = sorted(thermo, key=lambda x: (x[0], x[1]))
                # Save resonance isomers reordered by their thermo
                #spc.molecule = [item[2] for item in thermo]
            thermo0 = thermo[0][3] 
        else:
            # Saturated molecule, does not need HBI method
            thermo0 = thermodatabase.getThermoDataFromGroups(spc)
            
        # Make sure to calculate Cp0 and CpInf if it wasn't done already
        thermodatabase.findCp0andCpInf(spc, thermo0)
    return processThermoData(spc, thermo0, thermoClass)


def evaluator(aug_inchi):
    """
    Module-level function passed to workers.

    generates a thermodata object for the 
    identifier and stores it in the 
    thermo database.

    Next, thermo is generated of this species.

    """
    logging.debug("Evaluating augmented inchi %s ", aug_inchi)

    spc = Species(molecule=[Molecule().fromAugmentedInChI(aug_inchi)])
    thermo = generateThermoData(spc)
    transportData = generateTransportData(spc)

    return aug_inchi, thermo, transportData

def _callback(data):
    """Callback that retrieves data from future and updates database."""
     
    assert thermo_engine is not None

    identifier, thermodata, transportdata = data
    key = parse_inchi(identifier)

    logging.debug("Future of %s has called back.", identifier)
    if thermodata is not None and transportdata is not None:

        thermo_engine.update_db(key, (thermodata, transportdata))

    else:
        logging.error("Thermo for identifier %s is None...", identifier)
        raise Exception




class ThermoEngine(AbstractEngine):
    """docstring for ThermoEngine"""
    def __init__(self):
        super(ThermoEngine, self).__init__()
        self.liquid = LiquidModel()
    
    def parse(self, identifier):
        """
        Converts the identifier into a key suitable to serve
        as a database key.
        """
        return parse_inchi(identifier)

    def isValidID(self, identifier):
        """Check if the parameter string is a valid identifier"""
        return isValid(identifier)

    def call_back(self, data):
        _callback(data)

    def submit_future(self, identifier):
        data = evaluator(identifier)
        return data

class LiquidModel(object):
    """docstring for LiquidModel"""
    solventName = None
    solventData = None
    solventViscosity = None
    diffusionTemp = None

def generateTransportData(species):
    """
    Generate the transportData parameters for the species.
    """
    database = rmgpy.data.rmg.database
    transport_db = database.transport

    return transport_db.getTransportProperties(species)[0]
