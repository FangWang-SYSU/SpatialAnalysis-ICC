import numpy as np
from scipy.sparse import issparse, csr_matrix, csc_matrix
from glmpca.glmpca import glmpca
import warnings

def sparse_binomial_deviance(X, sz):
    
    # Calculate LP & L1P
    sz_expanded = np.repeat(sz, X.getnnz(axis=1))
    P = X.data / sz_expanded
    LP_data = np.log(P, where = P > 0)
    L1P_data = np.log1p(- P, where = P < 1)

    col_indices = X.indices
    term1 = np.bincount(col_indices, weights=X.data*(LP_data - L1P_data), minlength=X.shape[1])
    term2 = np.bincount(col_indices, weights=sz_expanded*L1P_data, minlength=X.shape[1])
    
    # Calculate ll_sat
    ll_sat = term1 + term2
    
    # 计算零模型
    feature_sums = X.sum(axis=0).A1
    sz_sum = sz.sum()
    p = feature_sums / sz_sum
    l1p = np.log1p(-p)
    ll_null = feature_sums * (np.log(p) - l1p) + sz_sum * l1p
    
    deviance = 2 * (ll_sat - ll_null)
    return deviance


def dense_binomial_deviance(X, sz):

    P = X / sz[:, None]
    L1P = np.log1p(-P)
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        log_P = np.log(P)
    
    ll_sat = np.nansum(X * (log_P - L1P) + sz[:, None] * L1P, axis=0)
    
    feature_sums = X.sum(axis=0)
    sz_sum = sz.sum()
    p = feature_sums / sz_sum
    l1p = np.log1p(-p)
    ll_null = feature_sums * (np.log(p) - l1p) + sz_sum * l1p
    
    deviance = 2 * (ll_sat - ll_null)
    return deviance


def sparse_poisson_deviance(X, sz):

    sz_expanded = np.repeat(sz, X.getnnz(axis=1))
    P = X.data / sz_expanded
    LP_data = np.log(P, where = P > 0)
    
    col_indices = X.indices
    ll_sat = np.bincount(col_indices, weights = X.data * LP_data, minlength = X.shape[1])
    
    feature_sums = X.sum(axis=0).A1
    sz_sum = sz.sum()
    ll_null = feature_sums * np.log(feature_sums / sz_sum)
    
    deviance = 2 * (ll_sat - ll_null)
    return deviance


def dense_poisson_deviance(X, sz):
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        log_ratio = np.log(X / sz[:, None])
    
    ll_sat = np.nansum(X * log_ratio, axis=0)
    feature_sums = X.sum(axis=0)
    sz_sum = sz.sum()
    ll_null = feature_sums * np.log(feature_sums / sz_sum)
    
    deviance = 2 * (ll_sat - ll_null)
    return deviance


def _compute_deviance_batch(m, fam='binomial', batch=None):

    if batch is None:
        return _compute_deviance(m, fam)
    else:
        if len(batch) != m.shape[0]:
            raise ValueError("The length of batch must match the number of cells.")

        levels = np.unique(batch)
        res = np.zeros((len(levels), m.shape[1]))
        for i in range(len(levels)):
            level = levels[i]
            res[i, :] = _compute_deviance(m[batch == level, :], fam)
        return res.sum(0).flatten()


def _compute_deviance(m, fam):

    sz = np.array(m.sum(axis=1)).flatten()
    if fam == 'poisson':
        lsz = np.log(sz)
        sz = np.exp(lsz - np.mean(lsz))  # Geometric mean normalization

    if issparse(m):
        if fam == 'binomial':
            deviance = sparse_binomial_deviance(m, sz)
        else:
            deviance = sparse_poisson_deviance(m, sz)
    else:
        if fam == 'binomial':
            deviance = dense_binomial_deviance(m, sz)
        else:
            deviance = dense_poisson_deviance(m, sz)
    
    deviance[np.isnan(deviance)] = 0
    return deviance


def _check_input(m):
    if not issparse(m) and isinstance(m, (csr_matrix, csc_matrix)):
        warnings.warn('X is not sparse matrix, convert to numpy array')
        m = m.toarray()
    elif issparse(m) and not isinstance(m, csr_matrix):
        m = csr_matrix(m)
        warnings.warn("The input matrix is sparse, converting to CSR matrix.")
    elif not issparse(m) and isinstance(m, np.ndarray):
        warnings.warn('X is not sparse matrix, convert to numpy array')
        int_m = m.astype(int)
        if np.array_equal(m, int_m):
            m = int_m
            print('X is closed to integer matrix, converting to integer matrix')
        else:
            warnings.warn("The input matrix is floating-point, integer values are recommended.")
    else:
        data = m.data
        if not np.issubdtype(m.dtype, np.integer):
            int_data = data.astype(int)
            if np.array_equal(data, int_data):
                m = csr_matrix(m.astype(int))
                print('X is closed to integer matrix, converting to integer matrix')
            else:
                warnings.warn("The input matrix is floating-point, integer values are recommended.")
    return m


def _devianceFeature(adata, fam = 'binomial', layer = None, batch = None, batch_key = None):
    # TODO: 
    # add the parameter for dense matrix option

    # Check input parameters
    if fam not in ['binomial', 'poisson']:
        raise ValueError("fam must be 'binomial' or 'poisson'")

    if layer is not None:
        m = adata.layers[layer].copy()
    else:
        m = adata.X.copy()
    
    m = _check_input(m)

    # Get batch information
    if batch is None and batch_key is not None:
        if batch_key not in adata.obs.columns:
            raise ValueError(f"Column '{batch_key}' not found in adata.obs")
        else:
            batch = adata.obs[batch_key].values
    elif batch is not None and batch_key is not None:
        raise ValueError("Cannot specify both batch and batch_key")

    # Compute the deviance
    deviance = _compute_deviance_batch(m, fam = fam, batch = batch)

    # Add the deviance to adata.var
    if 'deviance' in adata.uns.keys():
        if fam != adata.uns['deviance']['fam']:
            warnings.warn(f"The deviance family has changed from {adata.uns['deviance']['fam']} to {fam}.")
            adata.uns.pop('deviance')
    adata.uns['deviance'] = {'fam': fam, 'layer': layer, 'batch': batch, 'batch_key': batch_key}
    adata.var['deviance'] = deviance


def _run_glmpca(adata, fam = 'poi', n_components = None, n_genes = None, layer = None, penalty = 1, n_iter = 1000, eps = 1e-4, 
               verbose = False, **kwargs):

    if layer is not None:
        counts = adata.layers[layer].copy()
    else:
        counts = adata.X.copy()

    if type(counts).__module__.startswith('scipy.sparse'):
        counts = counts.toarray()

    if n_components is None:
        n_components = 10
        print("n_components not specified, using {} components".format(n_components))

    assert penalty > 0, "penalty must be positive"
    assert eps > 0, "eps must be a positive number"

    assert (isinstance(n_iter, int) and (n_iter > 0)), "n_iter must be an positive integer"
    assert (isinstance(n_components, int) and (n_components > 0)), "n_components must be an positive integer"

    gene_idx = adata.var['deviance'].values.argsort()[::-1]
    if n_genes is None:
        print("n_genes not specified, using all genes")
    else:
        print("n_genes specified, using top {} genes".format(n_genes))
        counts = counts[:, gene_idx[:n_genes]]
    
    return glmpca(counts.T, n_components, fam = fam, penalty = penalty, verbose = verbose, 
                  ctl = {'maxIter': n_iter, 'eps': eps, 'optimizeTheta': True}, **kwargs)


# TODO: 
# add a function to plot the results
# add a function to make DE analysis
# Integrated the function into A user-difined Class