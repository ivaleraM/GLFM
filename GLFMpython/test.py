from aux import plot_dim
import numpy as np

N = 50
D = 3
Kest = 2
maxR = 2
X = np.random.rand(D,N) + 3
X[1,:] = np.random.randint(1,5,N)
X[2,:] = np.random.randint(1,3,N) # 2 categories
B = np.zeros((D,Kest,maxR))
B[:,:,0]  = np.array([[3, 3],[3,1],[3,-1]]) # (D,Kest,maxR)
B[:,:,1]  = np.array([[0, 0],[0,0],[3,-1]]) # (D,Kest,maxR)
Theta = []
C = 'pnc'
d = 2
k = 0
s2Y = 1
s2u = 0
lab = ('A','B')
# Signature plot_dim(X,B,Theta,C,d,k,s2Y,s2u)
plot_dim(X,B,Theta,C,d,k,s2Y,s2u,labels=lab)
