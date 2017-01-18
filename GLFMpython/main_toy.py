import numpy as np
import sys

sys.path.append('../Ccode/wrapper_python/')

import GLFM

import pdb

Z = np.array([[1.0,0],[1,1],[1,1]])
X = np.array([[1.0, 1, -0.3, 1],[6, 2, 3.8, 23],[11, 3, 4.1, 4]])
C = 'PoGN'

print 'First, in Python'
print X
print Z

#print '\nX flags'
#print X.flags
#print '\n'
#print X.transpose().flags

#print '\nZ flags'
#print Z.flags
#print '\n'
#print Z.transpose().flags

#print '\nX flags'
#print X.flags
X2 = np.ascontiguousarray(X.transpose())
Z2 = np.ascontiguousarray(Z.transpose())
#print '\n'
#print X2.flags
#print X2

print '\nNow, inside C'
pdb.set_trace()
(Z_out,B_out,Theta_out) = GLFM.infer(X2,C,Z2)

print Z_out
print "\n"
print B_out
print "\n"
print Theta_out

print "SUCCESSFUL"

