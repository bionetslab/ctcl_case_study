import os
import scanpy as sc
import squidpy as sq
import numpy as np
import itertools as itt
import pandas as pd
from sklearn.metrics.pairwise import pairwise_distances
import networkx as nx
from scipy import sparse

if __name__ == '__main__':
    graph_types=['delaunay_graph', 'knn_graph', 'dist_thresh_graph']
    cnt=-1
    sq_mi_results=[]
    for filename in os.listdir('../../data'):
        filename = os.fsdecode(filename)
        if filename.endswith('.h5ad'):
            cnt+=1
            one_or_zero_flag=-1
            adata = sc.read_h5ad('../../data/'+filename)
            celltypes=list(sorted(set(adata.obs['celltype'])))
            nshape=len(list(adata.to_df().columns))
            mi_scores_per_graph_type=[]
            for graph_type in graph_types:
                if graph_type=='delaunay_graph':
                    sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)
                elif graph_type=='knn_graph':
                    samples = adata.obsm['spatial']
                    from sklearn.neighbors import NearestNeighbors
                    neigh = NearestNeighbors(n_neighbors=2)
                    neigh.fit(samples)
                    nearest_neighbors=[]
                    nearest_neighbor_distances=[]
                    for i in samples:
                        nearest_neighbor_distances.append(neigh.kneighbors([i])[0][0][1])
                        nearest_neighbors.append(neigh.kneighbors([i])[1][0][1])
                    G_nearest_neighbors=nx.Graph()
                    G_nearest_neighbor_distances=nx.Graph()
                    edges_nearest_neighbors=[(list(pd.DataFrame(nearest_neighbors).index)[i], list(pd.DataFrame(nearest_neighbors)[0])[i]) for i in range(0, len(list(pd.DataFrame(nearest_neighbors).index)))]
                    edges_nearest_neighbor_distances=[(list(pd.DataFrame(nearest_neighbor_distances).index)[i], list(pd.DataFrame(nearest_neighbor_distances)[0])[i]) for i in range(0, len(list(pd.DataFrame(nearest_neighbor_distances).index)))]
                    G_nearest_neighbors.add_edges_from(edges_nearest_neighbors)
                    weighted_edges_knn_graph=[(list(pd.DataFrame(nearest_neighbors).index)[i], list(pd.DataFrame(nearest_neighbors)[0])[i], list(pd.DataFrame(nearest_neighbor_distances)[0])[i]) for i in range(0, len(list(pd.DataFrame(nearest_neighbor_distances).index)))]
                    G_nearest_neighbor_distances.add_weighted_edges_from(weighted_edges_knn_graph)
                    A_nearest_neighbors = nx.to_numpy_array(G_nearest_neighbors)
                    A_nearest_neighbors=sparse.csr_matrix(A_nearest_neighbors)
                    A_nearest_neighbor_distances = nx.to_numpy_array(G_nearest_neighbor_distances, weight='weight')
                    A_nearest_neighbor_distances=sparse.csr_matrix(A_nearest_neighbor_distances)
                    adata.obsp['spatial_connectivities']=A_nearest_neighbors
                    adata.obsp['spatial_distances']=A_nearest_neighbor_distances
                elif graph_type=='dist_thresh_graph':
                    X = adata.obsm['spatial']
                    Y = adata.obsm['spatial']
                    pairwise_distances_matrix=pairwise_distances(X, Y, metric='sqeuclidean')
                    # # Pairwise distance matrix filtering for the computation of adjacency matrix:
                    # # i. Filtering by absolute value:
                    # # pairwise_distances_matrix_filtered = np.where(pairwise_distances_matrix > $a$, 0, pairwise_distances_matrix), where $a$: a real number
                    # # ii. Filtering by measures of central tendency (mean, median):
                    # # pairwise_distances_matrix_filtered = np.where(pairwise_distances_matrix > $a$, 0, pairwise_distances_matrix), where $a$: {np.mean(pairwise_distances_matrix), np.median(pairwise_distances_matrix)}
                    # # iii. Filtering by lowest $a$ percentile of distance values [i.e., each node pair ${n1, n2}$ in the adjacency matrix pairwise_distances_matrix with Euclidean distance $d(n1,n2)$ betwen them smaller than the lowest $a$ percentile of all pairwise distances, shall be retained, all others dropped]:
                    # # pairwise_distances_matrix_filtered = np.where(pairwise_distances_matrix > np.percentile(pairwise_distances_matrix, $a$), 0, pairwise_distances_matrix), where $a$: $R \in \[0,100\]$
                    A_nearest_neighbor_distance_matrix = np.where(pairwise_distances_matrix > np.percentile(pairwise_distances_matrix, 5), 0, pairwise_distances_matrix)
                    A_nearest_neighbor_distances = sparse.csr_matrix(A_nearest_neighbor_distance_matrix)
                    A_nearest_neighbor_matrix = np.where(A_nearest_neighbor_distance_matrix > 0, 1, 0)
                    A_nearest_neighbors = sparse.csr_matrix(A_nearest_neighbor_matrix)
                    adata.obsp['spatial_connectivities']=A_nearest_neighbors
                    adata.obsp['spatial_distances']=A_nearest_neighbor_distances
                sq.gr.spatial_autocorr( adata, mode="moran", genes=list(adata.to_df().columns), n_perms=100, n_jobs=1)
                moranI=adata.uns['moranI']
                moranI['protein']=list(moranI.index)
                n_shape=np.shape(moranI)[0]
                moranI['patient_id']=[adata.uns['patient_id'][0]]*nshape
                moranI['condition']=[adata.uns['Group'][0]]*nshape
                mi_scores_per_graph_type.append(moranI)
            mi_scores_per_graph_type=pd.concat(mi_scores_per_graph_type, axis=1)
    sq_mi_results.append(mi_scores_per_graph_type)
    sq_mi_results=pd.concat(sq_mi_results, axis=0)
    sq_mi_results.to_csv(os.path.join('../../results', 'squidpy_moranI_results.csv'))
    
    
    
    
    