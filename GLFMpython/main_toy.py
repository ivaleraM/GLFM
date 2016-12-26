import numpy as np
import gsl_run as gs

X = np.array([[1.0,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15]])
Z = np.array([[1.0,0],[1,1],[1,1]])
C = 'PPPPP'

print 'First, in Python'
print X
print Z

#print '\nX flags'
#print X.flags
#print '\n'
#print X.transpose().flags
#
#print '\nZ flags'
#print Z.flags
#print '\n'
#print Z.transpose().flags

#X2 = np.ascontiguousarray(X.transpose())

print '\nNow, inside C'
#gs.wrapper_IBPsampler(X,C,Z)
(Z,B,Theta) = gs.wrapper_IBPsampler(X.transpose(),C,Z.transpose())

print Z
print "\n"
print B
print "\n"
print Theta
print "SUCCESSFUL"

