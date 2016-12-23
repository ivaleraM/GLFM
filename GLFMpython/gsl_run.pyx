from cython_gsl cimport *
import cython
#from libc.stdlib cimport malloc, free
from cymem.cymem cimport Pool

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

cdef extern from "gsl/gsl_sf_exp.h":
    double gsl_sf_exp(const double x)

# declare the interface to the C code
#cdef extern void c_multiply (double* array, double value, int m, int n)
cdef extern from "stdio.h":
    int tolower(int c)

cdef extern from "InferenceFunctions.h":
    #void IBPsampler_toy(double missing, gsl_matrix* X, char *C, gsl_matrix* Z,\
    #                gsl_matrix **B)
    int IBPsampler_func (double missing, gsl_matrix *X, char *C, gsl_matrix *Z, gsl_matrix **B, gsl_vector **theta, int *R, double *w, int maxR, int bias, int N, int D, int K, double alpha, double s2B, double s2Y,int maxK,int Nsim);

@cython.boundscheck(False)
@cython.wraparound(False)
def wrapper_IBPsampler(np.ndarray[double, ndim=2, mode="c"] input not None,\
        Cin, np.ndarray[double, ndim=2, mode="c"] Zin not None, int bias=1,\
        double s2Y=1.0, double s2B=1.0, double alpha=1.0, int Nsim=50,\
        int maxK=10, double missing=-1):#\
    """
    """

    # cdef gsl_matrix * B = gsl_matrix_alloc(dim, m)
    cdef int N, D, K
    cdef gsl_matrix_view Xview, Zview
    cdef gsl_matrix* X
    cdef gsl_matrix* Zm
    cdef gsl_matrix** B

    N, D = input.shape[0], input.shape[1]
    K = Zin.shape[0] # TODO: Check input dims of Zin
    print 'N=%d, D=%d, K=%d\n' % (N, D, K)

    if len(Cin) != D:
        raise Exception('Size of C and X are not consistent!')
    #print 'size of C = 1*%d' % len(C)

    Xview = gsl_matrix_view_array(&input[0,0],N,D) #D,N)
    X = &Xview.matrix

    Zview = gsl_matrix_view_array(&Zin[0,0], K,N) # TODO: Check input dims
    Zm = &Zview.matrix # we need to allocate input matrix Z to [maxK*N] matrix
    cdef gsl_matrix* Z = gsl_matrix_calloc(maxK,N)
    for i in xrange(N):
        for k in xrange(K):
             gsl_matrix_set (Z, k, i,gsl_matrix_get (Zm, k, i))

    #B = <gsl_matrix **>malloc(sizeof(gsl_matrix*));
    cdef Pool mem = Pool()
    B = <gsl_matrix**>mem.alloc(D, sizeof(gsl_matrix*));
    cdef gsl_vector** theta = <gsl_vector**>mem.alloc(D, sizeof(gsl_vector*));

    ##...............BODY CODE.......................##

    C = ''
    for d in xrange(D):
        C += chr( tolower(ord(Cin[d])) ) # convert to lower case

    cdef np.ndarray[np.int32_t, ndim=1, mode="c"] R = np.empty(D,dtype=np.int32)
    cdef np.ndarray[double, ndim=1, mode="c"] w = np.empty(D)
    cdef np.ndarray[double, ndim=1] maxX = np.empty(D)
    cdef int maxR = 1
    #print 'Size of R[D]=(%d,%d)' % (R.shape[0],R.shape[1])
    #print 'Size of w[D]=(%d,%d)' % (w.shape[0],w.shape[1])
    #print 'Size of maxX[D]=(%d,%d)' % (maxX.shape[0],maxX.shape[1])
    #R[0] = 0; R[1] = 1; R[2] = 2; R[221] = 3
    #print 'R[0]=%d, R[1]=%d, R[2]=%d, R[3]=%d\n' % (R[0], R[1], R[2], R[221])

    cdef gsl_vector_view Xd_view
    for d in xrange(D):
        Xd_view = gsl_matrix_column(X,d)
        maxX[d] = gsl_vector_max(&Xd_view.vector)
        R[d] = 1
        w[d] = 1
        if C[d] == 'g':
            B[d] = gsl_matrix_alloc(maxK,1)
        elif C[d] == 'p':
            B[d] = gsl_matrix_alloc(maxK,1)
            w[d]=2/maxX[d]
        elif C[d] == 'n':
            B[d] = gsl_matrix_alloc(maxK,1)
            w[d] = 2/maxX[d]
        elif 'c':
            R[d] = <int>maxX[d]
            B[d] = gsl_matrix_alloc(maxK,R[d])
            if (R[d] > maxR):
                maxR = R[d]
        elif 'o':
            R[d] = <int>maxX[d];
            B[d] = gsl_matrix_alloc(maxK,1)
            theta[d] = gsl_vector_alloc(R[d])
            if (R[d] > maxR):
                maxR = R[d]


    #IBPsampler_toy(missing, X, C, Z, B) #&X[0,0]) #, value, m, n)
    cdef int Kest = IBPsampler_func(missing, X, C, Z, B, theta,\
            <int*> R.data, &w[0],\
            maxR, bias,  N, D, K, alpha, s2B, s2Y, maxK, Nsim);

    return None
