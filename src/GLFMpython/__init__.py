''' bnpy module __init__ file
'''
import os
import sys

# Configure PYTHONPATH before importing any bnpy modules
root = os.path.sep.join(os.path.abspath(__file__).split(os.path.sep)[:-2])
sys.path.append(os.path.join(root, 'datasets/'))
sys.path.append(os.path.join(root, 'third-party/'))
sys.path.append(os.path.join(root, 'third-party/anchorwordtopics/'))

import data
import suffstats
import util

import allocmodel
import obsmodel
from HModel import HModel

import ioutil
import init
import learnalg
import birthmove
import mergemove
import deletemove

import Run

run = Run.run
load_model = ioutil.ModelReader.load_model
save_model = ioutil.ModelWriter.save_model

__all__ = ['run', 'Run', 'learnalg', 'allocmodel', 'obsmodel', 'suffstats',
                   'HModel', 'init', 'util', 'ioutil']

## Configure save location
#hasWriteableOutdir = False
#if 'BNPYOUTDIR' in os.environ:
#    outdir = os.environ['BNPYOUTDIR']
#    if os.path.exists(outdir):
#        try:
#            testfilepath = os.path.join(
#                outdir, '.bnpy-test-write-permissions.txt')
#            with open(testfilepath, 'w') as f:
#                pass
#        except IOError:
#            sys.exit('BNPYOUTDIR not writeable: %s' % (outdir))
#        hasWriteableOutdir = True
#if not hasWriteableOutdir:
#    raise ValueError(
#        'Environment variable BNPYOUTDIR not specified.' +
#        ' Cannot save results to disk')
#
## Configure custom dataset directory
#if 'BNPYDATADIR' in os.environ:
#    if os.path.exists(os.environ['BNPYDATADIR']):
#        sys.path.append(os.environ['BNPYDATADIR'])
#    else:
#        print "Warning: Environment variable BNPYDATADIR not a valid directory"

# Optional viz package for plotting
try:
    from matplotlib import pylab
    import viz
    canPlot = True
    __all__.append('viz')
except ImportError:
    print "Error importing matplotlib. Plotting disabled."
    print "Fix by making sure this produces a figure window on your system"
    print " >>> from matplotlib import pylab; pylab.figure(); pylab.show();"
