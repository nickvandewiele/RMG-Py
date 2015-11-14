import os.path
import logging

from rmgpy.tools.loader import loadRMGPyJob


def loadReductionInput(reductionFile):
    """
    Load an reduction job from the input file located at `reductionFile`
    """

    target = None
    tolerance = -1

    full_path = os.path.abspath(os.path.expandvars(reductionFile))
    try:
        f = open(full_path)
    except IOError, e:
        logging.error('The input file "{0}" could not be opened.'.format(full_path))
        logging.info('Check that the file exists and that you have read access.')
        raise e

    logging.info('Reading input file "{0}"...'.format(full_path))
    
    global_context = { '__builtins__': None }
    local_context = {
        '__builtins__': None,
        'target': target,
        'tolerance': tolerance
    }

    try:
        exec f in global_context, local_context

        target = local_context['target']
        tolerance = local_context['tolerance']

    except (NameError, TypeError, SyntaxError), e:
        logging.error('The input file "{0}" was invalid:'.format(full_path))
        logging.exception(e)
        raise
    finally:
        f.close()

    assert target is not None
    assert tolerance != -1

    return target, tolerance

def load(rmg_inputFile, reductionFile, chemkinFile, speciesDict):
    
    rmg = loadRMGPyJob(rmg_inputFile, chemkinFile, speciesDict, generateImages=False)
    target, tolerance = loadReductionInput(reductionFile)

    return rmg, target, tolerance