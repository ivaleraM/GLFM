import numpy as np
import gsl_run as gs

X = np.array([[1.0,2,3,4,5],[6,7,8,9,10]])
Z = np.array([[1.0,0],[1,1]])
C = 'GPNCO'

print 'First, in Python'
print X
print Z.transpose()

print '\nNow, inside C'
gs.wrapper_IBPsampler(X,C,Z)

print "SUCCESSFUL"

