import numpy as np
import gsl_run as gs

Z = np.array([[1.0,0],[1,1],[1,1]])
X = np.array([[1.0, 1, -0.3, 1],[6, 2, 3.8, 23],[11, 3, 4.1, 4]])
C = 'PoGN'

#X = np.array([[1.0, 2, 3, 4, -0.3],[6, 7, 8, 9, 3.8],[11, 12, 13, 14, 4.1]])
#C = 'PPPPG' #'PPPPP'

#X = np.array([[1.0,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15]])
#C = 'PPPPG'

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

#print '\nX flags'
#print X.flags
X2 = np.ascontiguousarray(X.transpose())
Z2 = np.ascontiguousarray(Z.transpose())
#print '\n'
#print X2.flags
#print X2

print '\nNow, inside C'
#gs.wrapper_IBPsampler(X,C,Z)
(Z_out,B_out,Theta_out) = gs.wrapper_IBPsampler(X2,C,Z2)

print Z_out
print "\n"
print B_out
print "\n"
print Theta_out

print "SUCCESSFUL"

