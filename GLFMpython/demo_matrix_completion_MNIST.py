## import necessary packages
import numpy as np # library to work with numpy arrays and math operations
from random import sample
import sys
sys.path.append('../Ccode/wrapper_python/') # add relative path (location of GLFM library)
import GLFM

import pdb
import csv
import matplotlib.pyplot as plt
import matrix_completion as MC # import General Latent Feature Model Package

## load input data + add perc. of missings
print 'Loading data...'
# DB can be downloaded here: http://pjreddie.com/projects/mnist-in-csv/
file = '../databases/mnist_train.csv'
images = []
with open(file, 'r') as csv_file:
    count = 0
    for data in csv.reader(csv_file):
        count = count + 1
        # The first column is the label
        label = data[0]

        # The rest of columns are pixels
        pixels = np.array(data[1:], dtype='float64')
        images.append(pixels)

        # Make those columns into a array of 8-bits pixels
        # This array will be of 1D with length 784
        # The pixel intensity values are integers from 0 to 255
        pixels = np.array(pixels, dtype='uint8')

        if count > 1000:
            break

X = np.array(images).transpose() # D*N

N = 100
Xtrue = X[:,sample(xrange(X.shape[1]),N)] + 1.0 # add one, since 'n' type cannot start in zero
C = np.tile('g',(1,Xtrue.shape[0]))[0].tostring()
print 'Add missing values...'

perc_missing = 0.2
missing_val = -1
mask_missing = np.random.rand(Xtrue.shape[0],Xtrue.shape[1]) < perc_missing
Xmiss = Xtrue
Xmiss[mask_missing] = missing_val

## visualization of a random image
print 'Visualizing a single example...'
# Reshape the array into 28 x 28 array (2-dimensional array)
pixels = X[:,np.random.randint(0,X.shape[1])]
pixels = np.array(pixels, dtype='uint8')
pixels = pixels.reshape((28, 28))
# Plot
plt.imshow(pixels, cmap='gray',interpolation='none')
plt.ion() # interactive mode for plotting (script continues)
plt.show()

## Run algorithm (inference + estimation of missing values)
print 'Complete matrix...'
Kinit = 10
Z = np.ascontiguousarray( np.random.randint(0,2,size=(Kinit,N)).astype('float64') )
Xmiss = np.ascontiguousarray(Xmiss)

#pdb.set_trace()
#(Z_out,B_out,Theta_out) = GLFM.infer(Xmiss,C,Z)
#pdb.set_trace()
Xcompl = MC.complete_matrix(Xmiss, C, Niter=5, missing=missing_val) #, bias=0, s2Y=1, s2B=1, alpha=1, Niter=50, missing=-1)

