import numpy as np
import random
import sys
sys.path.append('../Ccode/wrapper_python/')

import scipy.io
import csv

class Data:

    def __init__(self, X, C, xlabel=[], ylabel=[], cat_labels=[]):
        """
        X: observation matrix, ´D x N´ array
        C: datatype, string of length ´1 x D´
        missing: value for missing data
        """
        self.X = X
        self.C = C
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.cat_labels = cat_labels

    @classmethod
    def load_from_csv(input_file):
        with open(input_file, 'rb') as csvfile:
            R = csv.reader(csvfile, delimiter=',')
            for row in R:
                # TODO: complete
                print ', '.join(row)
        data = Data(X,C,xlabel,ylabel,cat_labels)
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

        data = Data(X,C,xlabel,ylabel,cat_labels)
        return data

#    def preprocess(self):
#        """
#        Function to make matrix X suitable for Ccode inference routine
#        For count data, categorical and ordinal data, make sure that first
#        categories start at 1
#        For positive real data, make sure that the smallest value is bigger than 0
#
#        Inputs:
#           self: Data to be normalized
#
#        Output:
#           data: Data Object with normalized matrix X
#        """
#        (D,N) = self.X.shape
#        X2 = self.X.copy() # preprocess matrix
#        for d in xrange(D): # for each dimension
#            # get vector removing missing values in dimension d
#            if np.isnan(self.missing):
#                mask = np.where(not(np.isnan(self.X[d,:])))[0]
#            else:
#                mask = np.where(self.X[d,:] != self.missing)[0]
#            if mask.shape[0] == 0: # empty
#                continue
#            offset = min(self.[d,mask])
#            if self.C[d] == 'n' or self.C[d] == 'c' or self.C[d] == 'o':
#                X2[d,mask] = self.X[d,mask] - offset + 1
#            elif self.C[d] == 'p':
#                X2[d,mask] = self.X[d,mask] - offset + 10**-6
#            elif self.C[d] == 'g':
#                continue # TODO: Verify transformation for 'g'
#            else:
#                print 'Unkown datatype'
#        data = Data(X2,C,self.xlabel,self.ylabel,self.cat_labels,self.ylabel_long)
#        return data
