import sys
import os
#sys.path.append('../../src/GLFMpython/')
root = os.path.sep.join(os.path.abspath(__file__).split(os.path.sep)[:-3])
sys.path.append(os.path.join(root, 'src/GLFMpython/'))

import numpy as np
import GLFM

# --------------------------------------------------
# Function to verify good behavior of GLFM library
# --------------------------------------------------

print('\n\n# -------------------------------------------------')
print('# SCRIPT TO CHECK CALL TO GLFM LIBRARY')
print('# -------------------------------------------------\n')

data = dict()
data['X'] = np.array([[1.0, 1, -0.3, 1, 1],[6.3, 2, 3.8, 23, 1],[11, 3, 4.1, 4, 2]]) # (N,D)
data['C'] = 'poGNc'

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

