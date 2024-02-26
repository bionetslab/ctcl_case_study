import pandas as pd
import numpy as np

def _calculate_mean_column_value_per_category_(conditions, dataframe, category):
    df_mean_columns_values=[]
    for condition in conditions:
        df=dataframe[dataframe[category]==condition]
        mean_=df.mean(axis=0)
        mean_df=pd.DataFrame(mean_).T
        mean_df[category]=condition
        df_mean_columns_values.append(mean_df)
    df_mean_columns_values=pd.concat(df_mean_columns_values, axis=0)
    return

def _calculate_cell_results_by_sample_id_(celltypes, cell_results):
    cell_results_by_sample_id=[]
    for celltype in celltypes:
        # Filter cell_results data by celltype:
        cell_results_=cell_results[cell_results['cell_type']==celltype]
        
        # Prepare cell_results data by sample ids:
        sample_ids=list(set(cell_results_.sample_id))
        mean_local_entropy_3=[]
        mean_local_homophily_3=[]
        mean_egophily_3=[]
        condition=[]
        df_by_sample_ids=pd.DataFrame()
        df_by_sample_ids['sample_ids']=sample_ids
        for sampleId in sample_ids:
            sample_res=cell_results_[cell_results_['sample_id']==sampleId]
            mean_local_entropy_3.append(np.mean(list(sample_res['local_entropy_3'])))
            mean_local_homophily_3.append(np.mean(list(sample_res['local_homophily_3'])))
            mean_egophily_3.append(np.mean(list(sample_res['egophily_3'])))
            condition.append(list(sample_res.condition)[0])
        df_by_sample_ids['mean_local_entropy_3']=mean_local_entropy_3
        df_by_sample_ids['mean_local_homophily_3']=mean_local_entropy_3
        df_by_sample_ids['mean_egophily_3']=mean_egophily_3
        df_by_sample_ids['condition']=condition
        df_by_sample_ids['celltype']=celltype
        cell_results_by_sample_id.append(df_by_sample_ids)
    cell_results_by_sample_id=pd.concat(cell_results_by_sample_id, axis=0)
    return cell_results_by_sample_id

def _calculate_cell_results_by_patient_id_(celltypes, cell_results):
    cell_results_by_patient_id=[]
    for celltype in celltypes:
        # Filter cell_results data by celltype:
        cell_results_=cell_results[cell_results['cell_type']==celltype]
        
        # Prepare cell_results data by patient ids:
        patient_ids=list(set(cell_results_.patient_id))
        mean_local_entropy_3=[]
        mean_local_homophily_3=[]
        mean_egophily_3=[]
        condition=[]
        df_by_patient_ids=pd.DataFrame()
        df_by_patient_ids['patient_ids']=patient_ids.copy()
        patientIds_sampleIds_lookup={}
        for patientId in patient_ids:
            patient_res=cell_results_[cell_results_['patient_id']==patientId]
            patientIds_sampleIds_lookup[patientId]=list(set(patient_res['sample_id']))
            mean_local_entropy_3.append(np.mean(list(patient_res['local_entropy_3'])))
            mean_local_homophily_3.append(np.mean(list(patient_res['local_homophily_3'])))
            mean_egophily_3.append(np.mean(list(patient_res['egophily_3'])))
            condition.append(list(patient_res.condition)[0])
        df_by_patient_ids['mean_local_entropy_3']=mean_local_entropy_3
        df_by_patient_ids['mean_local_homophily_3']=mean_local_entropy_3
        df_by_patient_ids['mean_egophily_3']=mean_egophily_3
        df_by_patient_ids['condition']=condition
        df_by_patient_ids['celltype']=celltype
        cell_results_by_patient_id.append(df_by_patient_ids)
    cell_results_by_patient_id=pd.concat(cell_results_by_patient_id, axis=0)
    return cell_results_by_patient_id
