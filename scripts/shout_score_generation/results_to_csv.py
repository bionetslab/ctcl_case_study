import os
import scanpy as sc
import sys
import pandas as pd
sys.path.insert(0, '../../')


if __name__ == '__main__':
    graph_types=['delaunay_graph', 'knn_graph', 'dist_thresh_graph']
    sample_results = {
        'sample_id': [],
        'patient_id': [],
        'condition': [],
        'age': [],
        'sex': [],
        'number_of_cells': [],
        'SHouT_execution_time': [],
        'global_entropy': [],
        'global_homophily_delaunay_graph': [],
        'global_homophily_knn_graph': [],
        'global_homophily_dist_thresh_graph': []
    }
    score_types = ['local_entropy', 'egophily']
    cell_results = {f'{score_type}_{radius}': [] for score_type in score_types for radius in [1, 2, 3, 4, 5, 10, 20, 50, 100, 500]}
    score_types=[]
    for graph_type in graph_types:
        for radius in [1, 2, 3, 4, 5, 10, 20, 50, 100, 500]:
            cell_results[f'local_homophily_{radius}_{graph_type}']=[]
    score_types = ['local_homophily', 'local_entropy', 'egophily']
    cell_results['cell_type'] = []
    cell_results['condition'] = []
    cell_results['sample_id'] = []
    cell_results['patient_id'] = []

    for filename in os.listdir('../../results'):
        filename = os.fsdecode(filename)
        if filename.endswith('.h5ad'):
            adata = sc.read_h5ad('../../results/'+filename)
            sample_id = filename.split('.')[0]
            patient_id = adata.uns['patient_id'][0]
            condition = adata.uns['Group'][0]
            number_of_cells = len(adata.obs)
            # Update sample results
            sample_results['condition'].append(condition)
            sample_results['sample_id'].append(sample_id)
            sample_results['patient_id'].append(patient_id)
            sample_results['number_of_cells'].append(number_of_cells)
            sample_results['age'].append(int(adata.uns['Age'][0]))
            sample_results['sex'].append(adata.uns['Sex'][0])
            shout_results_list=['global_entropy']
            for graph_type in graph_types:
                shout_results_list.append(f'global_homophily_{graph_type}')
            shout_results_list.append('SHouT_execution_time')
            for shout_result in shout_results_list:
                sample_results[shout_result].append(adata.uns[shout_result])
            # Update cell results
            cell_results['condition'] += [condition] * number_of_cells
            cell_results['sample_id'] += [sample_id] * number_of_cells
            cell_results['patient_id'] += [patient_id] * number_of_cells
            cell_results['cell_type'] += list(adata.obs['celltype'])
            for score_type in score_types:
                for radius in [1, 2, 3, 4, 5, 10, 20, 50, 100, 500]:
                    if score_type=='local_homophily':
                        for graph_type in graph_types:
                            cell_results[f'{score_type}_{radius}_{graph_type}'] += list(adata.obs[f'{score_type}_{radius}_{graph_type}'])
                    else:
                        cell_results[f'{score_type}_{radius}'] += list(adata.obs[f'{score_type}_{radius}'])
    sample_results = pd.DataFrame(data=sample_results)
    cell_results = pd.DataFrame(data=cell_results)
    sample_results.to_csv(os.path.join('../../results', 'sample_results.csv'))
    cell_results.to_csv(os.path.join('../../results', 'cell_results.csv'))