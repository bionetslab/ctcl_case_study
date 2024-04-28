import os
import scanpy as sc
import sys
import pandas as pd
sys.path.insert(0, '../../')


if __name__ == '__main__':
    adata_all_samples_df=[]
    for filename in os.listdir('../data'):
        filename = os.fsdecode(filename)
        if filename.endswith('.h5ad'):
            adata = sc.read_h5ad('../data/'+filename)
            # Convert anndata object to pandas DataFrame:
            adata_df=adata.to_df()
            # Add .obs[] values as columns:
            adata_df=pd.concat([adata_df, adata.obs], axis=1)
            # Add .obsm[] values as columns:
            adata_df=pd.concat([adata_df, adata.obsm], axis=1)
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
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
            for shout_result in ['global_entropy', 'global_homophily', 'SHouT_execution_time']:
                sample_results[shout_result].append(adata.uns[shout_result])
            # Update cell results
            cell_results['condition'] += [condition] * number_of_cells
            cell_results['sample_id'] += [sample_id] * number_of_cells
            cell_results['patient_id'] += [patient_id] * number_of_cells
            cell_results['cell_type'] += list(adata.obs['celltype'])
            for score_type in score_types:
                for radius in [1, 2, 3, 4, 5]:
                    cell_results[f'{score_type}_{radius}'] += list(adata.obs[f'{score_type}_{radius}'])
    sample_results = pd.DataFrame(data=sample_results)
    cell_results = pd.DataFrame(data=cell_results)
    sample_results.to_csv(os.path.join('results', 'sample_results.csv'))
    cell_results.to_csv(os.path.join('results', 'cell_results.csv'))