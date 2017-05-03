import sys
sys.path.append('../../src/GLFMpython/')

import numpy as np
import GLFM
import pdb

# --------------------------------------------------
# Function to verify good behavior of GLFM library
# --------------------------------------------------

print('\n\n# -------------------------------------------------')
print('# SCRIPT TO CHECK CALL TO GLFM LIBRARY')
print('# -------------------------------------------------\n')

data = dict()
data['X'] = np.array([[1.0, 1, -0.3, 1],[6, 2, 3.8, 23],[11, 3, 4.1, 4]]) # (N,D)
data['C'] = 'goGN'

hidden = dict()
hidden['Z'] = np.array([[1.0,0],[1,1],[1,1]]) # dimensions (N,K)

print 'First, in Python'
print data['X']
print hidden['Z']

print '\nNow, inside C\n'
hidden = GLFM.infer(data,hidden)

print '\nBack to Python\n'
print hidden['Z']
print "\n"
print hidden['B']
print "\n"
print hidden['theta']
print hidden['mu']
print hidden['w']
print hidden['s2Y']

print('\n# -------------------')
print "# SUCCESSFUL"
print('# -------------------')

