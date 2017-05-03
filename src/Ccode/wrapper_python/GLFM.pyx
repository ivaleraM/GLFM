from cython_gsl cimport *
import cython
from cymem.cymem cimport Pool

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

import pdb

# declare the interface to the C code
#cdef extern void c_multiply (double* array, double value, int m, int n)
cdef extern from "stdio.h":
    int tolower(int c)

cdef extern from "../core/InferenceFunctions.h":
    int initialize_func (int N, int D, int maxK, double missing, gsl_matrix *X, char *C, gsl_matrix **B, gsl_vector **theta, int *R, double *f, double *mu,  double *w, double *s2Y)

cdef extern from "../core/InferenceFunctions.h":
    int IBPsampler_func (double missing, gsl_matrix *X, char *C, gsl_matrix *Z, gsl_matrix **B, gsl_vector **theta, int *R, double *f, double *mu, double *w, int maxR, int bias, int N, int D, int K, double alpha, double s2B, double *s2Y, double s2u, int maxK,int Nsim)
    # This is the C++ function to perform inference for the GLFM model
    # Inputs:
    #           missing: value of missings, cannot be nan, should be an integer
    #           X: observation matrix # TODO: Add matrix dimensions
    #           C: char vector to specify datatype of each dimension in X
    #           Z: binary matrix (feature activation matrix)
    #           B: feature matrix
    #           theta: auxiliary variables (necessary for ordinal)
    #           R: number of categories for each variable
    #           w: normalization weights in the transformation function # TODO: Make these weights input dependent
    #           maxR: maximum number of categories
    #           bias: how many columns should not be sampled
    #           N: number of observations
    #           D: number of dimensions
    #           K: number of latent features
    #           alpha: mass parameter of the IBP
    #           s2B: variance for feature values
    #           s2Y: variance for pseudo-observations
    #           s2u: auxiliary noise # TODO: Better explain
    #           maxK: maximum number of latent features (to allocate memory)
    #           Nsim: number of iterations (inside C++ code)
    # Outputs:
    #           Kest: number of inferred active features

@cython.boundscheck(False)
@cython.wraparound(False)
def infer(np.ndarray[double, ndim=2, mode="c"] Xin not None,\
        Cin not None, np.ndarray[double, ndim=2, mode="c"] Zin not None,\
        np.ndarray[double, ndim=1, mode="c"] Fin, int bias=0, double s2u=1.0,\
        double s2B=1.0, double alpha=1.0, int Nsim=100,\
        int maxK=50, double missing=-1, int verbose=0):#\
    """
    Function to call inference routine for GLFM model from Python code
    Inputs:
        Xin: observation matrix ( numpy array [D*N] )
        Cin: string array of length D
        Zin: latent feature binary matrix (numpy array [K*N] )
        Fin: vector of transform functions indicators

        *** (the following are optional parameters) ***
        bias: number of columns that should not be sampled in matrix Zin
        s2Y: variance for pseudo-observations Y
        s2u: auxiliary variance noise
        s2B: variance for feature values
        alpha: mass parameter for the IBP
        Nsim: number of iterations
        maxK: máximum number of latent features (for memory allocation)
        missing: value of missings (should be an integer or nan) # TODO: check
    Outputs:
        B_out: feature matrix: np.array of dimensions (D,Kest,maxR) where D is
               the number of dimensions, Kest is the number of inferred latent
               features, and maxR is the maximum number of categories
        Z_out: activation matrix: np.arry of dimensions (Kest,N) where Kest is
               the number of inferred latent features, and N = number of obs.
        theta_out: auxiliary variables for ordinal variables, ndarray of size
                   (D,maxR) where D = nr. of dimensions, maxR = max nr. of
                   categories
    """
    #print X_in.flags['C_CONTIGUOUS']
    #X_in = np.ascontiguousarray(X_in)
    #print X_in.flags['C_CONTIGUOUS']

    # cdef gsl_matrix * B = gsl_matrix_alloc(dim, m)
    cdef int N, D, K
    cdef gsl_matrix_view Xview, Zview
    cdef gsl_matrix* X
    cdef gsl_matrix* Zm

    N, D = Xin.shape[1], Xin.shape[0]
    K = Zin.shape[0]
    if verbose:
        print 'N=%d, D=%d, K=%d\n' % (N, D, K)
        print Xin
        print Zin

    ## transpose input matrices in order to be able to call inner C function
    #cdef np.ndarray[double, ndim=2, mode="c"] Xin = Xin_in.transpose()
    #cdef np.ndarray[double, ndim=2, mode="c"] Zin = Zin_in.transpose()

    if len(Cin) != D:
        raise Exception('Size of C and X are not consistent!')

    Zview = gsl_matrix_view_array(&Zin[0,0], K,N)
    Zm = &Zview.matrix # we need to allocate input matrix Z to [maxK*N] matrix
    cdef gsl_matrix* Z = gsl_matrix_calloc(maxK,N)
    for i in xrange(N):
        for k in xrange(K):
             gsl_matrix_set (Z, k, i,gsl_matrix_get (Zm, k, i))

    C = ''
    for d in xrange(D):
        C += chr( tolower(ord(Cin[d])) ) # convert to lower case

    ##...............BODY CODE.......................##

    Xview = gsl_matrix_view_array(&Xin[0,0],D,N)
    X = &Xview.matrix

    cdef Pool mem = Pool()
    cdef gsl_matrix** B = <gsl_matrix**>mem.alloc(D, sizeof(gsl_matrix*));
    cdef gsl_vector** theta = <gsl_vector**>mem.alloc(D, sizeof(gsl_vector*));

    cdef np.ndarray[double, ndim=1, mode="c"] w = np.empty(D)
    cdef np.ndarray[double, ndim=1, mode="c"] mu = np.empty(D)
    cdef np.ndarray[double, ndim=1, mode="c"] s2Y = np.empty(D)
    cdef np.ndarray[np.int32_t, ndim=1, mode="c"] R = np.empty(D,dtype=np.int32)
    print "In C++: transforming input data..."
    cdef int maxR = initialize_func(N, D, maxK, missing, X, C, B, theta,\
            <int*> R.data, &Fin[0], &mu[0], &w[0], &s2Y[0])
    print "done\n"
    #cdef int maxR = 1
    #cdef np.ndarray[double, ndim=1, mode="c"] Win = np.empty(D)

    if verbose:
        print "maxR = %d" % maxR

    ##...............Inference Function.......................##
    print '\nEntering C++: Running Inference Routine...\n'
    cdef int Kest = IBPsampler_func(missing, X, C, Z, B, theta,\
            <int*> R.data, &Fin[0], &mu[0], &w[0],\
            maxR, bias,  N, D, K, alpha, s2B, &s2Y[0], s2u, maxK, Nsim);
    print '\nBack to Python: OK\n'

    #print "w[0]=%.2f, w[1]=%.2f\n" % (float(w[0]), float(w[1]))

    ##...............Set Output Pointers.......................##
    cdef np.ndarray[double, ndim=2] Z_out = np.zeros((Kest,N))
    cdef np.ndarray[double, ndim=3] B_out = np.zeros((D,Kest,maxR))
    cdef np.ndarray[double, ndim=2] theta_out = np.zeros((D,maxR))
    #cdef np.ndarray[double, ndim=1] mu_out = np.zeros(D)
    #cdef np.ndarray[double, ndim=1] w_out = np.zeros(D)
    #cdef np.ndarray[double, ndim=1] s2Y_out = np.zeros(D)

    if verbose:
        print "Kest=%d, N=%d\n" % (Kest,N)

    # #Zview = gsl_matrix_submatrix(Z, 0, 0, Kest, N)
    #print 'N=%d, K=%d, Kest=%d' % (N,K,Kest)
    for i in xrange(N):
        for k in xrange(Kest):
            #print 'k=%d, i=%d' % (k,i)
            Z_out[k,i] = gsl_matrix_get(Z,k,i)
            #Z_out[k,i] = (&Zview.matrix)->data[k*N+i]
    if verbose:
        print "Z_out loaded"

    cdef gsl_matrix_view Bd_view
    cdef gsl_matrix* BT
    cdef int idx_tmp
    print "B_out[D,Kest,maxR] where D=%d, Kest=%d, maxR=%d" % (D,Kest,maxR)
    for d in xrange(D):
        #print 'd=%d, R[d]=%d' % (d,R[d])
        if (C[d] == 'o'):
            idx_tmp = 1
        else:
            idx_tmp = R[d]
        Bd_view =  gsl_matrix_submatrix(B[d], 0, 0, Kest, idx_tmp)
        BT = gsl_matrix_alloc(idx_tmp,Kest)
        gsl_matrix_transpose_memcpy (BT, &Bd_view.matrix)
        for k in xrange(Kest):
            #print 'k=%d' % k
            for i in xrange(idx_tmp):
                #print 'i=%d' % i
                B_out[d,k,i] = gsl_matrix_get(BT,i,k)
        gsl_matrix_free(BT)
    if verbose:
        print "B_out loaded"

    for d in xrange(D):
        for i in xrange(maxR):
            if (C[d]=='o' and i<(R[d]-1)):
                theta_out[d,i] = gsl_vector_get(theta[d],i)
    if verbose:
        print "theta_out loaded"

#    for d in xrange(D):
#        mu_out[d] = gsl_vector_get(mu,d)
#    if verbose:
#        print "mu_out loaded"
#
#    for d in xrange(D):
#        w_out[d] = gsl_vector_get(w,d)
#    if verbose:
#        print "w_out loaded"
#
#    for d in xrange(D):
#        s2Y_out[d] = gsl_vector_get(s2Y,d)
#    if verbose:
#        print "s2Y_out loaded"

    ##..... Free memory.....##
    for d in xrange(D):
        gsl_matrix_free(B[d])
        if (C[d] == 'o'):
            gsl_vector_free(theta[d])
    gsl_matrix_free(Z)

    return (Z_out,B_out,theta_out,mu, w, s2Y)

#    cdef gsl_vector_view Xd_view
#    for d in xrange(D):
#        Xd_view = gsl_matrix_row(X,d)
#        maxX[d] = gsl_vector_max(&Xd_view.vector)
#        if verbose:
#            print "maxX[%d] = %f\n" % (d,maxX[d])
#        R[d] = 1
#        #Win[d] = 1
#        if C[d] == 'g':
#            B[d] = gsl_matrix_alloc(maxK,1)
#        elif C[d] == 'p':
#            B[d] = gsl_matrix_alloc(maxK,1)
#            #Win[d]=2/maxX[d]
#        elif C[d] == 'n':
#            B[d] = gsl_matrix_alloc(maxK,1)
#            #Win[d] = 2/maxX[d]
#        elif C[d] == 'c':
#            R[d] = <int>maxX[d]
#            B[d] = gsl_matrix_alloc(maxK,R[d])
#            if (R[d] > maxR):
#                maxR = R[d]
#        elif C[d] == 'o':
#            R[d] = <int>maxX[d];
#            print 'Setting R[d]=%f' % R[d]
#            B[d] = gsl_matrix_alloc(maxK,1)
#            theta[d] = gsl_vector_alloc(R[d])
#            if (R[d] > maxR):
#                maxR = R[d]
