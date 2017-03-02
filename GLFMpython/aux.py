import numpy as np
import pdb
import mapping_functions as mf
import matplotlib.pyplot as plt

def preprocess(X,C,missing=-1):
    """
    Function to make matrix X suitable for Ccode inference routine
    For count data, categorical and ordinal data, make sure that first
    categories start at 1
    For positive real data, make sure that the smallest value is bigger than 0
    Inputs:
             X: observation matrix, D*N array
             C: datatype, string of length D
       missing: (optional) value for missing data
    Output:
       X2: standardize data matrix
    """
    (D,N) = X.shape
    X2 = X.copy() # preprocess matrix
    for d in xrange(D): # for each dimension
        # get vector removing missing values in dimension d
        if np.isnan(missing):
            mask = np.where(not(np.isnan(X[d,:])))[0]
        else:
            mask = np.where(X[d,:] != missing)[0]
        if mask.shape[0] == 0: # empty
            continue
        offset = min(X[d,mask])
        if C[d] == 'n' or C[d] == 'c' or C[d] == 'o':
            X2[d,mask] = X[d,mask] - offset + 1
        elif C[d] == 'p':
            X2[d,mask] = X[d,mask] - offset + 10**-6
        elif C[d] == 'g':
            continue
            # TODO
        else:
            print 'Unkown datatype'
    return X2

def plot_dim(X,B,Theta,C,d,k,s2Y,s2u,missing=-1,catlabel=[],xlabel=[]):
    """
    Function to plot an individual dimension of feature matrix B
    Inputs:
        X: observation matrix of dimensions (D,N)
        B: (D,Kest,maxR) ndarray
        C: datatype vector - str of length D
        d: dimension to plot
        k: feature to consider
    Output:
        void
    """
    plt.figure()
    plt.xlabel(xlabel)

    (D,Kest,maxR) = B.shape
    Xd = X[d,:]
    Cd = C[d]
    if k<0 or k>Kest:
        print('Error: k index should be bigger than o and smaller than Kest')
    if np.isnan(missing):
        mask = np.where(not(np.isnan(Xd)))[0]
    else:
        mask = np.where(Xd != missing)[0]
    if Cd == 'g':
        numP = 100 # number of points to plot
        xx = np.linspace( min(Xd[mask]), max(Xd[mask]), numP )
        Zn = np.zeros(Kest)
        Zn[k] = 1
        Bdv = B[d,:,0]
        pdf = mf.pdf_real(xx, Zn,Bdv,s2Y,s2u)
        plt.plot(xx,pdf)
        plt.ion()
        plt.show()

    elif Cd == 'p':
        numP = 100 # number of points to plot
        xx = np.linspace( min(Xd[mask]), max(Xd[mask]), numP )
        Zn = np.zeros(Kest)
        Zn[k] = 1
        Bdv = B[d,:,0]
        w = 2.0 / max(Xd[mask]) # TODO: put function handler
        pdf = mf.pdf_pos(xx,Zn,Bdv,w,s2Y,s2u,lambda x,w: mf.fpos_1(x,w), \
                lambda x,w: mf.dfpos_1(x, w));
        plt.plot(xx,pdf)
        plt.ion()
        plt.show()

    elif Cd == 'n':
        xx = np.arange( min(Xd[mask]), max(Xd[mask])+1)
        Zn = np.zeros(Kest)
        Zn[k] = 1
        Bdv = B[d,:,0]
        w = 2.0 / max(Xd[mask]) # TODO: put function handler
        pdf = mf.pdf_count(xx,Zn,Bdv,w,s2Y, lambda x,w: mf.fpos_1(x,w))
        plt.stem(xx,pdf)
        plt.ion()
        plt.show()

    elif Cd == 'c':
        R = len( np.unique(Xd[mask]) )
        Zn = np.zeros(Kest)
        Zn[k] = 1
        Bdv = np.squeeze(B[d,:,:]) # TODO: Review that size = [K*R]
        pdf = mf.pdf_cat(Zn,Bdv,s2u,R)
        bar_width = 0.35
        index = np.arange(len(pdf))
        plt.bar(index,pdf,bar_width)
        plt.xticks(index + bar_width / 2, catlabel, rotation='vertical')
        plt.ion()
        plt.show()

    elif Cd == 'o':
        a = 1
        # TODO
    else:
        print 'Unknown datatype'
    return
