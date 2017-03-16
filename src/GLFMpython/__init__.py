''' GLFM module __init__ file
'''
import os
import sys

# Configure PYTHONPATH before importing any bnpy modules
path = '/home/melanie/Documents/UC3M/Workspace/GIBP_Isabel/GLFM/'
sys.path.append(path+'src/Ccode/wrapper_python/')
root = os.path.sep.join(os.path.abspath(__file__).split(os.path.sep)[:-2])
print 'root=' + root
sys.path.append(os.path.join(root, 'Ccode/wrapper_python/'))
print os.path.join(root, 'Ccode/wrapper_python/')
print os.path.join(root, 'GLFMpython/')
#import data

__all__ = ['GLFM']
#__all__ = ['run', 'Run', 'learnalg', 'allocmodel', 'obsmodel', 'suffstats',
#                   'HModel', 'init', 'util', 'ioutil']

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

## Optional viz package for plotting
#try:
#    from matplotlib import pylab
#    import viz
#    canPlot = True
#    __all__.append('viz')
#except ImportError:
#    print "Error importing matplotlib. Plotting disabled."
#    print "Fix by making sure this produces a figure window on your system"
#    print " >>> from matplotlib import pylab; pylab.figure(); pylab.show();"
