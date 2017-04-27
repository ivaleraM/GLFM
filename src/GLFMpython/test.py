from GLFM import plot_dim
from GLFM import plot_dim_1feat
import numpy as np

N = 50
D = 3
Kest = 2
maxR = 2
X = np.random.rand(D,N) + 3
X[1,:] = np.random.randint(1,5,N)
X[2,:] = np.random.randint(1,3,N) # 2 categories
B = np.zeros((D,Kest,maxR))
B[:,:,0]  = np.array([[3, 2.5],[2,1],[2.4,-1]]) # (D,Kest,maxR)
B[:,:,1]  = np.array([[0, 0],[0,0],[3,-1]]) # (D,Kest,maxR)
Theta = []
C = 'gnc'
d = 0
k = 0
s2Y = 1
s2u = 0
lab = ('A','B')
# Signature plot_dim(X,B,Theta,C,d,k,s2Y,s2u)
#plot_dim_1feat(X,B,Theta,C,d,k,s2Y,s2u,catlabel=lab)

Zp = np.eye(Kest)
plot_dim(X,B,Theta,C,d,Zp,s2Y,s2u,catlabel=lab,xlabel='bla')
