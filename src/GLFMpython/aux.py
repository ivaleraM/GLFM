import numpy as np
import mapping_functions as mf
import matplotlib.pyplot as plt

#import cPickle
#input_file = '../databases/dataExploration/csv_xls/prostate.csv'
#with open(input_file) as f:
#    reader = csv.reader(f, delimiter=',')
#    for row in reader:
#        print row
#        pdb.set_trace()
#    [Y,genetic_ids, clinical_ids,vocab] = cPickle.load(f)

def preprocess(X,C,missing=-1):
    """
    Function to make matrix X suitable for Ccode inference routine
    For categorical and ordinal data, make sure that first
    categories start at 1
    For real-valued, positive and count data, make transformation
    For positive real data, make sure that the smallest value is bigger than 0
    For count data, smallest value should be 1
    Inputs:
             X: observation matrix, D*N array
             C: datatype, string of length D
       missing: (optional) value for missing data
    Output:
       X2: standardize data matrix
    """
    (D,N) = X.shape
    X2 = X.copy() # preprocess matrix
    suffStats = []
    for d in xrange(D): # for each dimension
        # get vector removing missing values in dimension d
        if np.isnan(missing):
            mask = np.where(not(np.isnan(X[d,:])))[0]
        else:
            mask = np.where(X[d,:] != missing)[0]
        if mask.shape[0] == 0: # empty
            continue
        if C[d] == 'g':
            mu = np.mean(X2[d,mask])
        elif C[d] == 'p':
            mu = min(X2[d,mask]) - 10**-10
        elif C[d] == 'n':
            mu = min(X2[d,mask]) - 1
        elif (C[d] == 'c') or (C[d] == 'o'):
            mu = min(X2[d,mask]) - 1
        else:
            print 'Unkown datatype'
        if (C[d] == 'g') or (C[d] == 'p') or (C[d] == 'n'):
            X2[d,mask] = 2.0 * (X[d,mask] - mu) / max(X[d,mask]-mu)
            suffStats.append(np.array([ -mu, 2.0/max(xs-mu) ]))
        else:
            X2[d,mask] = X[d,mask] - mu
            suffStats.append(np.array( [-mu] ))
    return (X2,suffStats)

def plot_dim_1feat(X,B,Theta,C,d,k,s2Y,s2u,missing=-1,catlabel=[],xlabel=[]):
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

    elif Cd == 'n':
        xx = np.arange( min(Xd[mask]), max(Xd[mask])+1)
        Zn = np.zeros(Kest)
        Zn[k] = 1
        Bdv = B[d,:,0]
        w = 2.0 / max(Xd[mask]) # TODO: put function handler
        pdf = mf.pdf_count(xx,Zn,Bdv,w,s2Y, lambda x,w: mf.fpos_1(x,w))
        plt.stem(xx,pdf)

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

    elif Cd == 'o':
        a = 1
        # TODO
    else:
        print 'Unknown datatype'
    plt.ion()
    plt.show()
    plt.pause(0.0001)
    return

def plot_dim(X,B,Theta,C,d,Zp,s2Y,s2u,missing=-1,catlabel=[],xlabel=[]):
    """
    Function to plot an individual dimension of feature matrix B
    Inputs:
        X: observation matrix of dimensions (D,N)
        B: (D,Kest,maxR) ndarray
        C: datatype vector - str of length D
        d: dimension to plot
    Output:
        void
    """
    if (Zp.shape[1] != B.shape[1]):
        print 'Error: Sizes of Zp and B are inconsistent'

    colors = ['r','b','g','m','g']
    plt.figure()       # create new figure
    plt.xlabel(xlabel) # add x legend
    #print xlabel

    (D,Kest,maxR) = B.shape
    Xd = X[d,:]
    Cd = C[d]
    # only consider values in dimension d which are not missing
    if np.isnan(missing):
        mask = np.where(not(np.isnan(Xd)))[0]
    else:
        mask = np.where(Xd != missing)[0]
    (numPatterns,Kest) = Zp.shape
    if Cd == 'g':
        numP = 100 # number of points to plot
        xx = np.linspace( min(Xd[mask]), max(Xd[mask]), numP )
        Bdv = B[d,:,0]
        for p in xrange(numPatterns):
            Zn = np.squeeze(Zp[p,:]) # TODO: Verify dimensions
            pdf = mf.pdf_real(xx, Zn,Bdv,s2Y,s2u)
            plt.plot(xx,pdf,label=str(Zn))

    elif Cd == 'p':
        numP = 100 # number of points to plot
        xx = np.linspace( min(Xd[mask]), max(Xd[mask]), numP )
        Bdv = B[d,:,0]
        w = 2.0 / max(Xd[mask]) # TODO: put function handler
        for p in xrange(numPatterns):
            Zn = np.squeeze(Zp[p,:]) # TODO: Verify dimensions
            pdf = mf.pdf_pos(xx,Zn,Bdv,w,s2Y,s2u,lambda x,w: mf.fpos_1(x,w), \
                lambda x,w: mf.dfpos_1(x, w))
            plt.plot(xx,pdf,colors[p],label=str(Zn))

    elif Cd == 'n':
        xx = np.arange( min(Xd[mask]), max(Xd[mask])+1)
        Bdv = B[d,:,0]
        w = 2.0 / max(Xd[mask]) # TODO: put function handler
        for p in xrange(numPatterns):
            Zn = np.squeeze(Zp[p,:]) # TODO: Verify dimensions
            pdf = mf.pdf_count(xx,Zn,Bdv,w,s2Y, lambda x,w: mf.fpos_1(x,w))
            plt.stem(xx,pdf,colors[p], label=str(Zn))

    elif Cd == 'c':
        R = len( np.unique(Xd[mask]) )
        Bdv = np.squeeze(B[d,:,:]) # TODO: Review that size = [K*R]
        bar_width = 0.6/numPatterns
        for p in xrange(numPatterns):
            Zn = np.squeeze(Zp[p,:]) # TODO: Verify dimensions
            pdf = mf.pdf_cat(Zn,Bdv,s2u,R)
            index = np.arange(len(pdf))
            plt.bar(index+p*bar_width,pdf,width=bar_width,color=colors[p]) #, label=str(Zn))
        plt.xticks(index + bar_width / 2, catlabel, rotation='vertical')
#ax.bar(x-0.2, y,width=0.2,color='b',align='center')
#ax.bar(x, z,width=0.2,color='g',align='center')


    elif Cd == 'o':
        print 'This category is currently under development'
        a = 1
        # TODO
    else:
        print 'Unknown datatype'
    plt.legend()
    plt.ion()
    plt.show()
    plt.pause(0.0001)

    return
