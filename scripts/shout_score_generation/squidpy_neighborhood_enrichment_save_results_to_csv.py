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
    sq_ne_results=[]
    for filename in os.listdir('../../data'):
        filename = os.fsdecode(filename)
        if filename.endswith('.h5ad'):
            cnt+=1
            one_or_zero_flag=-1
            adata = sc.read_h5ad('../../data/'+filename)
            celltypes=list(sorted(set(adata.obs['celltype'])))
            cluster_pairs_list=list(itt.combinations(celltypes, 2))
            nshape_interactions=len(cluster_pairs_list)
            ne_scores_per_graph_type=[]
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
                nhood_enrichment=sq.gr.nhood_enrichment(adata, cluster_key='celltype', copy=True, backend="multiprocessing", show_progress_bar=False)
                nhood_enrichment_zscore=nhood_enrichment[0]
                nhood_enrichment_count=nhood_enrichment[1]
                n_shape=np.shape(nhood_enrichment_zscore)[0]
                upper_zscore_matrix = np.triu(np.array(nhood_enrichment_zscore), 1)
                zscore_list=[]
                for j in range(1, n_shape):
                    for k in range(j, n_shape):
                        zscore_list.append(upper_zscore_matrix[j,k])
                upper_count_matrix = np.triu(np.array(nhood_enrichment_count), 1)
                count_list=[]
                for j in range(1, n_shape):
                    for k in range(j, n_shape):
                        count_list.append(upper_count_matrix[j,k])
                df=pd.DataFrame(cluster_pairs_list, columns =['celltype_1', 'celltype_2'])
                df['sq_ne_zscore']=zscore_list
                df['sq_ne_count']=count_list
                df['patient_id']=[adata.uns['patient_id'][0]]*nshape_interactions
                df['condition']=[adata.uns['Group'][0]]*nshape_interactions
                ne_scores_per_graph_type.append(df)
            ne_scores_per_graph_type=pd.concat(ne_scores_per_graph_type, axis=1)
            sq_ne_results.append(ne_scores_per_graph_type)
    sq_ne_results=pd.concat(sq_ne_results, axis=0)
    sq_ne_results.to_csv(os.path.join('../../results', 'squidpy_nhoodEnrichment_results.csv'))
    
    
    
    
    