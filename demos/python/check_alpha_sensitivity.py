import numpy as np # import numpy matrix for calculus with matrices
import sys
sys.path.append('../../src/GLFMpython/')
import GLFM        # import General Latent Feature Model Library
import matplotlib.pyplot as plt # import plotting library
import time        # import time to be able to measure iteration speed

import pdb
from IPython import embed

# import libraries for I/O of data
import cPickle as pickle
import cPickle

import pdb
N_ITER_TRAINING = 10000
N_ITER_BETWEEN_SAMPLES = 30


def run_prostate(
        n_samples=None,
        alpha=None,
        seed=42,
        n_iter_training=N_ITER_TRAINING,
        n_iter_between_samples=N_ITER_BETWEEN_SAMPLES,
        maxK=50,
        ):

    with open('../../datasets/prostate.pk','rb') as f:
        data = cPickle.load(f)

    # ---------------------------------------------
    # 2. INITIALIZATION FOR GLFM ALGORITHM
    # ---------------------------------------------
    print('\n 2. INITIALIZATION\n')

    print('\tSetting optional parameters for the GLFM model...')

    [N, D] = data['X'].shape

    params = dict()

    # pre-transform a subset of variables
    idx_transform = [ D-1 ] # we transform the last dimension
    params['t'] = [None] * D
    params['t_1'] = [None] * D
    params['dt_1'] = [None] * D
    params['ext_dataType'] = [None] * D
    for rr in xrange(len(idx_transform)):
        r = idx_transform[rr]
        params['t_1'][r] = lambda x: np.log(x + 1) # transformation to apply to raw data
        params['t'][r] = lambda y: np.exp(y) - 1   # inverse transform to recover raw data
        params['dt_1'][r] = lambda x: 1/(x+1)      # derivative of inverse transform
        params['ext_dataType'][r] = 'p'    # change type of data due to transformation

    params['s2u'] = .005    # Auxiliary variance
    params['s2B'] = 1       # Variance of the Gaussian prior of the weigting matrices B
    params['maxK'] = maxK   # maximum number of latent features for memory allocation
    params['bias'] = 1      # 1 = fix first feature to be active for all patients 

    print('\tInitializing Z...')

    np.random.seed(seed)

    Kinit = 2   # initial number of latent features
    prob = 0.2  # probability of feature activation in matrix Z
    hidden = dict()
    if params['bias']:
        initZ = np.concatenate((np.ones((N,1)),(np.random.rand(N,Kinit-1) < 0.2)*1.0),axis=1)
    else:
        initZ = (np.random.rand(N,Kinit) < prob) * 1.0

    params['Niter'] = n_iter_training

    params['alpha'] = alpha     # Concentration parameter of the IBP
    hidden['Z'] = initZ
    # ---------------------------------------------
    # 3. RUN INFERENCE FOR GLFM ALGORITHM
    # ---------------------------------------------
    print('\n 3. INFERENCE\n')

    print('\tInfering latent features...')
    hidden = GLFM.infer(data,hidden,params=params)

    # ---------------------------------------------
    # 4. COLLECT SAMPLES
    # ---------------------------------------------
    params['Niter'] = n_iter_between_samples

    samples_Z = []
    for sample in range(n_samples):
        hidden = GLFM.infer(data,hidden,params=params)
        samples_Z.append(hidden['Z'])

    return samples_Z

def run_counties(
        n_samples=None,
        alpha=None,
        seed=42,
        n_iter_training=N_ITER_TRAINING,
        n_iter_between_samples=N_ITER_BETWEEN_SAMPLES,
        maxK=20,
        ):

    # ---------------------------------------------
    # 1. LOAD DATA TO BE EXPLORED
    # ---------------------------------------------
    print('\n 1. LOAD DATABASE TO EXPLORE\n')

    with open('../../datasets/counties.pk','rb') as f:
        data = cPickle.load(f)

    # ---------------------------------------------
    # 2. INITIALIZATION FOR GLFM ALGORITHM
    # ---------------------------------------------
    print('\n 2. INITIALIZATION\n')

    print('\tSetting optional parameters for the GLFM model...')

    params = dict()
    params['s2B'] = 1      # noise variance for feature values
    params['s2u'] = 0.005  # auxiliary noise

    [N, D] = data['X'].shape

    # pre-transform a subset of variables
    idx_transform = [ 1, 4, 5, 10] # dimensions to be transformed
    params['t'] = [None] * D
    params['t_1'] = [None] * D
    params['dt_1'] = [None] * D
    params['ext_dataType'] = [None] * D
    for rr in xrange(len(idx_transform)):
        r = idx_transform[rr]
        params['t_1'][r] = lambda x: np.log(x + 1) # transformation to apply to raw data
        params['t'][r] = lambda y: np.exp(y) - 1   # inverse transform to recover raw data
        params['dt_1'][r] = lambda x: 1/(x+1)      # derivative of inverse transform
        params['ext_dataType'][r] = 'p'    # change type of data due to transformation
    # dimension 'White' need an inversion too
    r = 9
    params['t_1'][r] = lambda x: np.log((100-x)+1)  # transformation to apply to raw data
    params['t'][r] = lambda y: - np.exp(y) + 101    # inverse transform to recover raw data
    params['dt_1'][r] = lambda x: -1/(101 - x)      # derivative of inverse transform
    params['ext_dataType'][r] = 'p'                 # change type of data due to transformation

    params['maxK'] = maxK     # maximum number of latent features for memory allocation
    params['bias'] = 1      # 1 = fix first feature to be active for all patients 


    print('\tInitializing Z...')

    np.random.seed(seed)

    Kinit = 2   # initial number of latent features
    prob = 0.2  # probability of feature activation in matrix Z
    hidden = dict()
    if params['bias']:
        initZ = np.concatenate((np.ones((N,1)),(np.random.rand(N,Kinit-1) < 0.2)*1.0),axis=1)
    else:
        initZ = (np.random.rand(N,Kinit) < prob) * 1.0

    params['Niter'] = n_iter_training

    params['alpha'] = alpha     # Concentration parameter of the IBP
    hidden['Z'] = initZ
    # ---------------------------------------------
    # 3. RUN INFERENCE FOR GLFM ALGORITHM
    # ---------------------------------------------
    print('\n 3. INFERENCE\n')

    print('\tInfering latent features...')
    hidden = GLFM.infer(data,hidden,params=params)

    # ---------------------------------------------
    # 4. COLLECT SAMPLES
    # ---------------------------------------------
    params['Niter'] = n_iter_between_samples

    samples_Z = []
    for sample in range(n_samples):
        hidden = GLFM.infer(data,hidden,params=params)
        samples_Z.append(hidden['Z'])

    return samples_Z


''' computes list of Z matrices of length n_samples; then save as npy file '''
def compute_GLFM_posterior_samples(
        dataset=None,
        output_path=None,
        alpha=1.0,
        seed=42,
        n_samples=500,
        n_iter_training=N_ITER_TRAINING,
        n_iter_between_samples=N_ITER_BETWEEN_SAMPLES,
        ):

    if dataset == 'counties':
        samples_Z = run_counties(
                n_samples=n_samples,
                alpha=alpha,
                seed=seed,
                n_iter_training=n_iter_training,
                n_iter_between_samples=n_iter_between_samples,
                )

    elif dataset == 'prostate':
        samples_Z = run_prostate(
                n_samples=n_samples,
                alpha=alpha,
                seed=seed,
                n_iter_training=n_iter_training,
                n_iter_between_samples=n_iter_between_samples,
                )
    else:
        raise ValueError('Unknown dataset')

    # save structure
    with open(output_path,'w') as f:
        pickle.dump(samples_Z,f)
    print('DONE')

if __name__ == '__main__':

    compute_GLFM_posterior_samples(
        dataset='counties',
        output_path='test.pkl',
        alpha=1.0,
        seed=42,
        n_samples=3,
        n_iter_training=10,
        n_iter_between_samples=2,
        )
