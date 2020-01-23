import numpy as np
import pickle

import matplotlib.pyplot as plt

import pdb

seed_vec = [0]#, 1, 2]
alpha_vec = [0.1, 1.0, 10.0, 100.0]
dataset = 'counties'
xlims=[0,20]

#seed_vec = [0, 1]
#dataset = 'prostate'
#xlims=[0,51]

perc_vec = [0.0, 0.001, 0.025, 0.05, 0.075, 0.01, 0.1]

#@cache
def load_Zs(output_path):
    with open(output_path,'r') as f:
        hidden_Z = pickle.load(f)
    return hidden_Z

for perc in perc_vec:

    f, ax = plt.subplots(4,1,figsize=(6,4), sharex=True)

    for aa, alpha in enumerate(alpha_vec):
        K_vec_all = []
        for seed in seed_vec:

            output_path = 'RESULTS/db={}-seed={}-alpha={}.pkl'.format(
                    dataset,seed,alpha)
            hidden_Z = load_Zs(output_path)
            print('LOADED {}'.format(output_path))

            nk = [np.sum(x,axis=0) for x in hidden_Z]
            N = hidden_Z[0].shape[0]

            # filtering
            K_vec = [np.sum(x > perc*N) for x in nk]
            K_vec_all += K_vec

        ax[aa].hist(K_vec, bins=np.linspace(0.5,50.5,51),
                density=True, label=('alpha = %.1f' % alpha))
        ax[aa].legend()
        ax[aa].set_xlim(xlims)
        if aa == 3:
            ax[aa].set_xlabel('Number of active features')

    plt.savefig('RESULTS/db={}-perc={}.png'.format(
        dataset,perc),bbox_inches='tight')

    #pdb.set_trace()
