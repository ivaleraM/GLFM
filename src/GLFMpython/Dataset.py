import numpy as np
import random
import sys
sys.path.append('../Ccode/wrapper_python/')

import scipy.io
import csv

class Dataset:

    def __init__(self, X, C, xlabel=[], ylabel=[], cat_labels=[],\
            ylabel_long=[],missing=-1):
        """
        X: observation matrix, D*N array
        C: datatype, string of length D
        missing: value for missing data
        """
        self.X = X
        self.C = C
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.cat_labels = cat_labels
        self.missing = missing

    @classmethod
    def load_from_csv(input_file):
        data = Dataset(X,C,xlabel,ylabel,cat_labels,ylabel_long)
        return data

    @classmethod
    def load_from_mat(input_file):
        tmp = scipy.io.loadmat(input_file)
        data = tmp['data'][0,0] # data is a dictionary with the following keys
        (N,D) = data['X'].shape
        X = data['X'].transpose() #  ndarray of dimensions D * N
        C = str(data['C'][0])
        # dealing with missing data: replace np.nan by -1
        (xx,yy) = np.where(np.isnan(X)) # find positions where X is nan (i.e. missing data)
        for r in xrange(len(xx)):
            X[xx[r],yy[r]] = -1

        data = Dataset(X,C,xlabel,ylabel,cat_labels,ylabel_long)
        return data

    def preprocess(self):
        """
        Function to make matrix X suitable for Ccode inference routine
        For count data, categorical and ordinal data, make sure that first
        categories start at 1
        For positive real data, make sure that the smallest value is bigger than 0

        Inputs:
           self: Dataset to be normalized

        Output:
           data: Dataset Object with normalized matrix X
        """
        (D,N) = self.X.shape
        X2 = self.X.copy() # preprocess matrix
        for d in xrange(D): # for each dimension
            # get vector removing missing values in dimension d
            if np.isnan(self.missing):
                mask = np.where(not(np.isnan(self.X[d,:])))[0]
            else:
                mask = np.where(self.X[d,:] != self.missing)[0]
            if mask.shape[0] == 0: # empty
                continue
            offset = min(self.[d,mask])
            if self.C[d] == 'n' or self.C[d] == 'c' or self.C[d] == 'o':
                X2[d,mask] = self.X[d,mask] - offset + 1
            elif self.C[d] == 'p':
                X2[d,mask] = self.X[d,mask] - offset + 10**-6
            elif self.C[d] == 'g':
                continue # TODO: Verify transformation for 'g'
            else:
                print 'Unkown datatype'
        data = Dataset(X2,C,self.xlabel,self.ylabel,self.cat_labels,self.ylabel_long)
        return data
