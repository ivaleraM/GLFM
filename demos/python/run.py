import argparse

from check_alpha_sensitivity import compute_GLFM_posterior_samples

N_SAMPLES=500
N_ITER_TRAINING = 10000
N_ITER_BETWEEN_SAMPLES = 30

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset',type=str,default='counties')

    parser.add_argument('--alpha',type=float,default=1.0)
    parser.add_argument('--seed',type=int,default=1)
    parser.add_argument('--n_samples',type=int,default=N_SAMPLES)
    parser.add_argument('--n_iter_training',type=int,default=N_ITER_TRAINING)
    parser.add_argument('--n_iter_between_samples',type=int,
            default=N_ITER_BETWEEN_SAMPLES)

    arg_dict = vars(parser.parse_args())
    locals().update(arg_dict)

    output_path = 'RESULTS/db={}-seed={}-alpha={}.pkl'.format(dataset,seed,alpha)

    compute_GLFM_posterior_samples(
        dataset=dataset,
        output_path=output_path,
        alpha=alpha,
        seed=seed,
        n_samples=n_samples,
        n_iter_training=n_iter_training,
        n_iter_between_samples=n_iter_between_samples,
        )
