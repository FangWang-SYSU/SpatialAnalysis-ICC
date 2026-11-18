from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Literal

import numpy as np
from scipy.sparse import issparse, csr_matrix, csc_matrix

if TYPE_CHECKING:
    from anndata import AnnData
    from numpy import ndarray

def sparse_binomial_deviance(X: csr_matrix, sz: ndarray):
    
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


def dense_binomial_deviance(X: ndarray, sz: ndarray):

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


def sparse_poisson_deviance(X: csr_matrix, sz: ndarray):

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


def dense_poisson_deviance(X: ndarray, sz: ndarray):
    
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
    if m.shape[1] > 8192:
        m = m.toarray()
    return m


def _devianceFeature(
        adata: AnnData, 
        fam: Literal['binomial', 'poisson'] = 'binomial', 
        layer: str | None = None, 
        batch: ndarray | None = None, 
        batch_key = None
    ):
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
    if batch is None:
        deviance = _compute_deviance(m, fam=fam)
    else:
        batch = np.array(batch)
        if len(batch) != m.shape[0]:
            raise ValueError("The length of batch must match the number of cells.")
        deviance = _compute_deviance_batch(m, fam=fam, batch=batch)

    #deviance = _compute_deviance_batch(m, fam = fam, batch = np.array(batch))

    # Add the deviance to adata.var
    if 'deviance' in adata.uns.keys():
        if fam != adata.uns['deviance']['fam']:
            warnings.warn(f"The deviance family has changed from {adata.uns['deviance']['fam']} to {fam}.")
            adata.uns.pop('deviance')
    adata.uns['deviance'] = {'fam': fam, 'layer': layer, 'batch': np.array(batch), 'batch_key': batch_key}
    adata.var['deviance'] = deviance


def deviance(adata, fam = 'binomial', layer = None, batch = None, batch_key = None):
    """
    Calculate the deviation for the AnnData object and add the result to AnnData
    
    Parameters
    ----------
    adata : AnnData
        Input AnnData object, assuming X or the specified layer is the count matrix (cells x genes)
    fam : str, optional
        Family of distributions, can be 'binomial' or 'poisson', default is 'binomial'
    layer : str, optional
        Use the specified layer instead of adata.X, default is None (use adata.X)
    batch : array-like, optional
        Batch information, default is None
    batch_key : str, optional
        Column name in adata.obs for specifying batch information, default is None
    
    Returns
    -------
    adata : AnnData
        Input AnnData object with the deviance added to adata.var
    """

    # Compute the deviance
    _devianceFeature(adata, fam, layer, batch, batch_key)

def highly_deviance_genes(
        adata: AnnData,
        n_genes: int = 5000,
    ):
    pass

# TODO: 
# add a function to plot the results
# add a function to make DE analysis
# Integrated the function into A user-difined Class