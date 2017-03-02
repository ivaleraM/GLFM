import numpy as np
import sys

sys.path.append('../Ccode/wrapper_python/')

import GLFM

import pdb

Z = np.array([[1.0,0],[1,1],[1,1]])
X = np.array([[1.0, 1, -0.3, 1],[6, 2, 3.8, 23],[11, 3, 4.1, 4]])
C = 'PoGN'
W = 2.0 / np.max(X,0)

print 'First, in Python'
print X
print Z

X2 = np.ascontiguousarray(X.transpose())
Z2 = np.ascontiguousarray(Z.transpose())

print '\nNow, inside C\n'
#pdb.set_trace()
(Z_out,B_out,Theta_out) = GLFM.infer(X2,C,Z2,W)

print '\nBack to Python\n'
print Z_out
print "\n"
print B_out
print "\n"
print Theta_out

print "SUCCESSFUL"

