from __future__ import annotations

import functools
import itertools
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal, Dict

import numpy as np
import pandas as pd
import scipy.sparse as sp
from glmpca.glmpca import glmpca as gpca

# Only run when type checking
if TYPE_CHECKING:
    from anndata import AnnData
    from numpy import ndarray


def retry_on_failure(max_retry: int = 3):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            attempts = 0
            last_exception = None
            
            for attempt in itertools.count(1):
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    last_exception = e
                    attempts += 1
                    print(f"Attempt {attempt}: {func.__name__} failed with {e}. Retrying...")
                    if attempts >= max_retry:
                        raise last_exception from None
        return wrapper
    return decorator


def encode(covar: ndarray, 
           encode_method: Literal['onehot', 'label'] = None):
    from sklearn.preprocessing import OneHotEncoder, LabelEncoder
    if not encode_method:
        encode_method = 'onehot'
    if encode_method == 'onehot':
        # One-hot encoding for categorical data
        encoder = OneHotEncoder(sparse_output=False)
        encoded_covar = encoder.fit_transform(covar.reshape(-1, 1))
    elif encode_method == 'label':
        # Label encoding for categorical data
        encoder = LabelEncoder()
        encoded_covar = encoder.fit_transform(covar).reshape(-1, 1)
    else:
        raise ValueError(f"Unsupported encode_method: {encode_method}")
    return encoded_covar


class _DataInput:

    adata: AnnData
    layer: str

    def __init__(
            self, 
            adata: AnnData,
            layer: str | None = None
        ) -> None:
        
        self._adata = adata
        
        if layer:
            Y = adata.layers[layer]
        else:
            Y = adata.X
        if type(Y).__module__.startswith('scipy.sparse'):
            self.Y = Y.toarray()
        else:
            self.Y = Y

    def get_batch_obs(
            self, 
            obs_covars: ndarray | None = None, 
            obs_keys: str | None = None, 
            obsm_keys: ndarray | None = None,
            encode_type: Dict[str, Literal['onehot', 'label']] = {}
        ):
        """
        Get the batch of observations for the GLM-PCA model.

        Parameters:
        -----------
        obs_covars : np.ndarray, optional
            Predefined observation covariates. If None, extract from adata.
        obs_keys : str or list, optional
            Keys for the observation covariates in adata.obs.
        obsm_keys : str or list, optional
            Keys for the observation covariates in adata.obsm.
        encode_type : dict, optional
            Encoding method for categorical data. Default is 'onehot' for all covariates.

        Attributes:
            X : np.ndarray (N x P)
                The batch of observations.
        """
        from scipy.sparse import issparse

        assert isinstance(obs_covars, (np.ndarray, type(None))), "obs_covars must be a numpy array."
        assert isinstance(obs_keys, (str, list, tuple, type(None))), "obs_keys must be a string or list."
        assert isinstance(obsm_keys, (str, list, tuple, type(None))), "obsm_keys must be a string or list."

        adata = self._adata
        covar_list = []

        if obs_covars is not None:
            if isinstance(obs_covars, np.ndarray):
                covar_list.append(obs_covars)
            else:
                raise ValueError("obs_covars must be a numpy array.")
            
        if obs_keys is not None:
            if isinstance(obs_keys, str):
                obs_keys = [obs_keys]
            for key in obs_keys:
                if key in adata.obs.columns:
                    covar = adata.obs[key].values
                    # Check if the data is categorical or string type
                    if covar.dtype.kind in ['O', 'S', 'U']:
                        encoded_covar = encode(covar, encode_type.get(key, 'onehot'))
                        covar_list.append(encoded_covar)
                    else:
                        covar_list.append(covar.reshape(-1, 1))
                else:
                    raise KeyError(f"Key {key} not found in adata.obs.")
        
        if obsm_keys is not None:
            if isinstance(obsm_keys, str):
                obsm_keys = [obsm_keys]
            for key in obsm_keys:
                if key in adata.obsm_keys():
                    covar = adata.obsm[key]
                    if issparse(covar):
                        covar = covar.toarray()
                    covar_list.append(covar)
                else:
                    raise KeyError(f"Key {key} not found in adata.obsm.")

        # Concatenate all covariates into one matrix
        if covar_list:
            X = np.hstack(covar_list)
        else:
            X = None

        self.X = X

    def get_batch_var(
            self, 
            var_covars: ndarray | None = None, 
            var_keys: str | None = None, 
            varm_keys: ndarray | None = None,
            encode_type: Dict[str, Literal['onehot', 'label']] = {}
        ):
        """
        Get the batch of observations for the GLM-PCA model.

        Parameters:
        -----------
        adata : AnnData
            The AnnData object containing the data.
        obs_covars : np.ndarray, optional
            Predefined observation covariates. If None, extract from adata.
        obs_keys : str or list, optional
            Keys for the observation covariates in adata.obs.
        obsm_keys : str or list, optional
            Keys for the observation covariates in adata.obsm.
        encode_type : dict, optional
            Encoding method for categorical data. Default is 'onehot' for all covariates.

        Attributes:
            Z : np.ndarray (GxP)
                The batch of observations.
        """
        from scipy.sparse import issparse

        assert isinstance(var_covars, (np.ndarray, type(None))), "obs_covars must be a numpy array."
        assert isinstance(var_keys, (str, list, tuple, type(None))), "obs_keys must be a string, list, or None."
        assert isinstance(varm_keys, (np.ndarray, tuple, type(None))), "obsm_keys must be a numpy array or None."

        adata = self._adata
        covar_list = []

        if var_covars is not None:
            if isinstance(var_covars, np.ndarray):
                covar_list.append(var_covars)
            else:
                raise ValueError("obs_covars must be a numpy array.")
            
        if var_keys is not None:
            if isinstance(var_keys, str):
                var_keys = [var_keys]
            for key in var_keys:
                if key in adata.var.columns:
                    covar = adata.var[key].values
                    # Check if the data is categorical or string type
                    if covar.dtype.kind in ['O', 'S', 'U']:
                        encoded_covar = encode(covar, encode_type.get(key, 'onehot'))
                        covar_list.append(encoded_covar)
                    else:
                        covar_list.append(covar.reshape(-1, 1))
                else:
                    raise KeyError(f"Key {key} not found in adata.var.")
        
        if varm_keys is not None:
            if isinstance(varm_keys, str):
                varm_keys = [varm_keys]
            for key in varm_keys:
                if key in adata.varm_keys():
                    covar = adata.varm[key]
                    if isinstance(covar, (sp.csr_matrix, sp.csc_matrix)):
                        covar = covar.toarray()
                    covar_list.append(covar)
                else:
                    raise KeyError(f"Key {key} not found in adata.varm.")

        # Concatenate all covariates into one matrix
        if covar_list:
            Z = np.hstack(covar_list)
        else:
            Z = None

        self.Z = Z

@dataclass
class _Params:

    fam: str = 'poi'
    npcs: int = 30
    penalty: int | float = 1
    nb_theta: int | float = 100
    n_iter: int = 100
    eps: int| float = 1e-4
    verbose: bool = False

    def __post_init__(self):
        assert isinstance(self.fam, str), "fam must be a string"
        assert self.fam in ['poi', 'binomial', 'multinomial', 'bernoulli'], "fam must be one of 'poi', 'binomial', 'multinomial', 'bernoulli'."
        assert isinstance(self.npcs, int) and self.npcs > 0, "npcs must be a positive integer"
        assert isinstance(self.penalty, (int, float)) and self.penalty >= 0, "penalty must be a non-negative number"
        assert isinstance(self.nb_theta, (int, float)) and self.nb_theta > 0, "nb_theta must be a positive number"
        assert isinstance(self.n_iter, int) and self.n_iter > 0, "n_iter must be a positive integer"
        assert isinstance(self.eps, (int, float)) and self.eps > 0, "eps must be a positive number"
        assert isinstance(self.verbose, bool), "verbose must be a boolean"

    def get_options(self):
        return {
            'L': self.npcs,
            'fam': self.fam,
            'ctl': {
                "maxIter": self.n_iter,
                "eps": self.eps,
                "optimizeTheta": True
            },
            'penalty': self.penalty,
            'verbose': self.verbose,
            'nb_theta': self.nb_theta
        }


def _glmpca(
        # Anndata object
        adata: AnnData, 
        layer: str = None, 
        obs_covars: ndarray | None = None, 
        obs_keys: str | None = None, 
        obsm_keys: ndarray | None = None,
        var_covars: ndarray | None = None, 
        var_keys: str | None = None, 
        varm_keys: ndarray | None = None,
        encode_type: Dict[str, Literal['onehot', 'label']] = {},
        deviance_key: str = 'deviance',
        # GLMPCA parameters
        fam: str = 'poi', 
        npcs: int = 30, 
        top_n: int = None,
        nb_theta = 100, 
        penalty = 1, 
        n_iter = 100, 
        eps = 1e-4, 
        verbose = False):

    # get input data
    data = _DataInput(adata, layer)
    data.get_batch_var(var_covars, var_keys, varm_keys, encode_type)
    data.get_batch_obs(obs_covars, obs_keys, obsm_keys, encode_type)
    Y = data.Y
    X = data.X
    Z = data.Z

    # get parameters
    params = _Params(fam, npcs, penalty, nb_theta, n_iter, eps, verbose)

    print("Calculate {} components".format(params.npcs))

    gene_idx = adata.var[deviance_key].values.argsort()[::-1]
    if top_n is None:
        print("n_genes not specified, using all genes")
    else:
        print("n_genes specified, using top {} genes".format(top_n))
        Y = Y[:, gene_idx[:top_n]]
        if Z is not None:
            Z = Z[gene_idx[:top_n], :]
    
    glmpca_result = gpca(Y.T, X = X, Z = Z, **params.get_options()) 
    config = params.get_options()
    config.update(
        {
            'dev': glmpca_result['dev'],
            'top_n_genes': top_n,
        }
    )
    return glmpca_result

def glmpca(
        # Anndata object
        adata: AnnData, 
        layer: str = None, 
        obs_covars: ndarray | None = None, 
        obs_keys: str | None = None, 
        obsm_keys: ndarray | None = None,
        var_covars: ndarray | None = None, 
        var_keys: str | None = None, 
        varm_keys: ndarray | None = None,
        encode_type: Dict[str, Literal['onehot', 'label']] = {},
        deviance_key: str = 'deviance',
        # GLMPCA parameters
        fam: str = 'poi', 
        npcs: int = 30, 
        top_n: int = None,
        nb_theta = 100, 
        penalty = 1, 
        n_iter = 100, 
        eps = 1e-4, 
        verbose = False, 
        n_retry = 3):
    """
    Run GLM-PCA on AnnData object with retry mechanism.

    Parameters:
    -----------
    adata : AnnData
        Input AnnData object.
    layer : str, optional
        Use the specified layer instead of adata.X, default is None (use adata.X).
    inplace : bool, optional
        Whether to modify the input AnnData object, default is False.
    copy : bool, optional
        Whether to return a copy of the input AnnData object, default is False.
    obs_covars : array-like, optional
        Observational covariates, default is None.
    obs_keys : str or array-like, optional
        Column names in adata.obs for specifying observational covariates, default is None.
    obsm_keys : array-like, optional
        Column names in adata.obsm for specifying observational covariates, default is None.
    var_covars : array-like, optional
        Variational covariates, default is None.
    var_keys : str or array-like, optional
        Column names in adata.var for specifying variational covariates, default is None.
    varm_keys : array-like, optional
        Column names in adata.varm for specifying variational covariates, default is None.
    encode_type : dict, optional
        Encoding type for categorical variables, default is {}.
    deviance_key : str, optional
        Column name in adata.var for specifying the deviance, default is 'deviance'.
    key_add : str, optional
        Key to add the result to adata.uns, default is 'glmpca'.
    fam : str, optional
        Family of GLM-PCA, default is 'poi'.
    npcs : int, optional
        Number of principal components to use, default is 30.
    top_n : int, optional
        Number of genes to use, default is None.
    nb_theta : int, optional
        Number of theta to use, default is 100.
    penalty : float, optional
        Penalty parameter, default is 1.
    n_iter : int, optional
        Number of iterations, default is 1000.
    eps : float, optional
        Epsilon value, default is 1e-4.
    verbose : bool, optional
        Whether to print verbose information, default is False.
    n_retry : int, optional
        Number of retries, default is 3.

    """
    @retry_on_failure(max_retry = n_retry)
    def _run_glmpca():
        return _glmpca(
            adata,
            layer = layer,
            obs_covars = obs_covars, 
            obs_keys = obs_keys, 
            obsm_keys = obsm_keys,
            var_covars = var_covars, 
            var_keys = var_keys, 
            varm_keys = varm_keys,
            encode_type = encode_type,
            deviance_key = deviance_key,
            # GLMPCA parameters
            fam = fam, 
            npcs = npcs, 
            top_n = top_n,
            nb_theta = nb_theta, 
            penalty = penalty, 
            n_iter = n_iter, 
            eps = eps, 
            verbose = verbose
        )
    glmpca_res = _run_glmpca()
    return glmpca_res