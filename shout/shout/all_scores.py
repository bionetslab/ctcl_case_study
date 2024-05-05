from .entropy import *
from .homophily import *
from .egophily import *
from .utility import *
from scipy.sparse.csgraph import shortest_path


def all_scores(adata, cluster_key, radii, normalize=True, num_cell_types=None, copy=False):
    """Saves as dict, or returns as part of input anndata object, all heterogeneity scores calculated using the SHouT package.

    Parameters
    ----------
    adata : anndata.AnnData (Mandatory parameter)
        Annotated data object containing spatial omics data with spatial coordinates stored in `adata`.obsm['spatial'].
        For more info, go to https://anndata.readthedocs.io/en/latest/.
    
    cluster_key : str (Mandatory parameter)
        adata.obs[cluster_key] contains key where clustering/ cell type annotations are stored.
    
    radii : list[int] (Mandatory parameter)
        List of n-hop neighbors over which local heterogeneity scores are to be calculated.
    
    normalize : bool (Optional parameter | default True)
        If $normalize = True$, normalize by number of cell types in adata.obs[cluster_key] when calculating entropy scores.
    
    num_cell_types : int (Optional parameter | default None)
        If $num_cell_types=None$, num_cell_types is the number of cell types present in adata.obs[cluster_key].
    
    copy : bool (Optional parameter | default False)
        $copy = True$ returns all scores as a dict, $copy = False$ saves all scores as part of the input anndata object "adata".
    
    Returns
    -------
    scores : dict [str, float/list[float]]
        Returns value if $copy = True$, saves to input anndata object "adata" if $copy = False$.

    """
    adj_matrix = get_spatial_graph(adata)
    shortest_path_distances = shortest_path(adj_matrix)
    adj_matrix_homophilic = get_homophilic_edges(adata, cluster_key, adj_matrix)
    scores = dict()
    scores['global_entropy'] = global_entropy(adata, cluster_key, normalize=normalize, num_cell_types=num_cell_types, copy=copy)
    scores['global_homophily'] = global_homophily(adata, cluster_key, copy=copy,
                                                  adj_matrix=adj_matrix, adj_matrix_homophilic=adj_matrix_homophilic)
    for radius in radii:
        extended_neighborhoods = get_extended_neighborhoods(shortest_path_distances, radius)
        scores[f'local_entropy_{radius}'] = local_entropy(adata, cluster_key, radius, normalize=normalize, num_cell_types=num_cell_types, 
                                                          copy=copy, shortest_path_distances=shortest_path_distances,
                                                          extended_neighborhoods=extended_neighborhoods)
        scores[f'local_homophily_{radius}'] = local_homophily(adata, cluster_key, radius, copy=copy,
                                                              adj_matrix=adj_matrix,
                                                              adj_matrix_homophilic=adj_matrix_homophilic,
                                                              shortest_path_distances=shortest_path_distances,
                                                              extended_neighborhoods=extended_neighborhoods)
        scores[f'egophily_{radius}'] = egophily(adata, cluster_key, radius, copy=copy,
                                                shortest_path_distances=shortest_path_distances,
                                                extended_neighborhoods=extended_neighborhoods)
    if copy:
        return scores
