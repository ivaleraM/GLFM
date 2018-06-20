import numpy as np
import mapping_functions as mf
from scipy.stats import norm

def compute_log_likelihood(
        Xtrue_ND,
        C,
        hidden,
        params,
        ):
    """
    This function calculates the log-lik of
    the held-out data in the GLFM. The function uses the trained
    parameters,
    to calculate the PDF at the points of missing values, using the true values.
    Inputs:
        - Xtrue_ND: N*D array of true data. The function computes the log-lik
                   values for all no nan values in X_true
        - C: D      string with data types
        - hidden:   inferred latent variables from the GLFM
        - params:   hyperparameters from the GLFM
    Output:
        - log-lik array
    """
    assert (hidden['Z'].shape[1] == hidden['B'].shape[1]),\
            "Incongruent sizes between Z and B"
    D = Xtrue_ND.shape[1]

    lik = np.zeros(X_true_ND.shape)

    # Apply external transformation if any
    Xtrue_transformed_ND = Xtrue_ND.copy()
    for d in xrange(D):
        if (params['t'][d] != None): # if there is an external transformation
            # change type of dimension d by external data type
            C = C[:d]+params['ext_dataType'][d]+C[d+1:]
            mask = (~np.nan(Xtrue_ND[:,d]))
            Xtrue_transformed_ND[mask,d] = params['t_1'][d](Xtrue_ND[mask,d])

    # Find coordinates of all missing values
    lik_cords = np.argwhere(~np.isnan(Xtrue_ND)).tolist()
    lik_cords = [tuple(l) for l in lik_cords]

    unique_maps = [None] * D
    for d in xrange(D):
        # Ensure that tnd is mapped to the correct low label
        if C[d] == 'c' or C[d] == 'o':
            unique_vals = np.unique(Xtrue_ND[:, d])
            # Categories have to start at 1
            unique_map = dict(zip(unique_vals, range(1, len(unique_vals) + 1)))
            unique_maps[d] = (unique_map)

    # Loop over all missing values, and calculate the log lik.
    for cord in lik_cords:
        n = cord[0]  # Get the row ID
        d = cord[1]  # Get the column ID
        tnd = Xtrue_ND[cord]  # Find the true value at the held out coordinate

        # ENSURE that tnd is mapped to the correct low label
        if C[d] == 'c' or C[d] == 'o':
            # Proper re-map
            tnd = unique_maps[d][tnd]

        if C[d] == 'g':
            # Real
            lik[cord] = mf.pdf_g(tnd,
                                            hidden['Z'][n, :],
                                            hidden['B'][d, :, 0],
                                            hidden['mu'][d],
                                            hidden['w'][d],
                                            hidden['s2Y'][d],
                                            params['s2u'])
        elif C[d] == 'p':
            # Positive real
            lik[cord] = mf.pdf_p(tnd,
                                            hidden['Z'][n, :],
                                            hidden['B'][d, :, 0],
                                            hidden['mu'][d],
                                            hidden['w'][d],
                                            hidden['s2Y'][d],
                                            params['s2u'])
        elif C[d] == 'n':
            # Count
            lik[cord] = mf.pdf_n(tnd,
                                            hidden['Z'][n, :],
                                            hidden['B'][d, :, 0],
                                            hidden['mu'][d],
                                            hidden['w'][d],
                                            hidden['s2Y'][d])
        elif C[d] == 'c':
            # Categorical
            pdf = mf.pdf_c(hidden['Z'][n, :],
                           hidden['B'][d, :, range(int(hidden['R'][d].shape[0]))],
                           hidden['s2Y'][d])
            lik[cord] = pdf[int(tnd - 1)]  # -1 since the max category is indexed by max_val -1

        elif C[d] == 'o':
            # Ordinal
            pdf = mf.pdf_o_single(int(tnd - 1), hidden['Z'][n, :],
                           hidden['B'][d, :, 0],
                           hidden['theta'][d, range(int(hidden['R'][d].shape[0] - 1))],
                           hidden['s2Y'][d])
            lik[cord] = pdf
            #log_lik[cord] = np.log(pdf[int(tnd - 1)])
        else:
            raise ValueError('Unknown data type')

    assert (np.isnan(lik).sum() == 0), "Some values are nan!"

    # transform lik pdf_y (pseudo-obs) into lik pdf_x
    for d in xrange(D):
        if params.has_key('t'):
            if (params['t'][d] != None): # we have used a special transform beforehand
                mask = (~np.nan(Xtrue_ND[:,d]))
                lik[mask,d] = lik[mask,d] * np.abs(
                    params['dt_1'][d](Xtrue_ND[mask,d]) )

    # finally, apply log transform
    for d in xrange(D):
        lik[mask,d] = np.log(lik[mask,d])
    return lik
