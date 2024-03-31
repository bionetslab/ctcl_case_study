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
    return df_mean_columns_values

def _calculate_cell_results_by_sample_id_(celltypes, cell_results):
    cell_results_by_sample_id=[]
    for celltype in celltypes:
        # Filter cell_results data by celltype:
        cell_results_=cell_results[cell_results['cell_type']==celltype]
        
        # Prepare cell_results data by sample ids:
        sample_ids=list(set(cell_results_.sample_id))
        mean_local_entropy_r=[]
        mean_local_homophily_r=[]
        mean_egophily_r=[]
        condition=[]
        df_by_sample_ids=pd.DataFrame()
        df_by_sample_ids['sample_ids']=sample_ids
        for sampleId in sample_ids:
            sample_res=cell_results_[cell_results_['sample_id']==sampleId]
            mean_local_entropy_r.append(np.mean(list(sample_res['local_entropy_5'])))
            mean_local_homophily_r.append(np.mean(list(sample_res['local_homophily_5'])))
            mean_egophily_r.append(np.mean(list(sample_res['egophily_5'])))
            condition.append(list(sample_res.condition)[0])
        df_by_sample_ids['mean_local_entropy_5']=mean_local_entropy_r
        df_by_sample_ids['mean_local_homophily_5']=mean_local_entropy_r
        df_by_sample_ids['mean_egophily_5']=mean_egophily_r
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
        mean_local_entropy_5=[]
        mean_local_homophily_5=[]
        mean_egophily_5=[]
        condition=[]
        df_by_patient_ids=pd.DataFrame()
        df_by_patient_ids['patient_ids']=patient_ids
        patientIds_sampleIds_lookup={}
        for patientId in patient_ids:
            patient_res=cell_results_[cell_results_['patient_id']==patientId]
            patientIds_sampleIds_lookup[patientId]=list(set(patient_res['sample_id']))
            mean_local_entropy_5.append(np.mean(list(patient_res['local_entropy_5'])))
            mean_local_homophily_5.append(np.mean(list(patient_res['local_homophily_5'])))
            mean_egophily_5.append(np.mean(list(patient_res['egophily_5'])))
            condition.append(list(patient_res.condition)[0])
        df_by_patient_ids['mean_local_entropy_5']=mean_local_entropy_5
        df_by_patient_ids['mean_local_homophily_5']=mean_local_homophily_5
        df_by_patient_ids['mean_egophily_5']=mean_egophily_5
        df_by_patient_ids['condition']=condition
        df_by_patient_ids['celltype']=celltype
        cell_results_by_patient_id.append(df_by_patient_ids)
    cell_results_by_patient_id=pd.concat(cell_results_by_patient_id, axis=0)
    return cell_results_by_patient_id

def _calculate_centrality_scores_by_patient_id_(celltypes, cell_results):
    cell_results_by_patient_id=[]
    for celltype in celltypes:
        # Filter cell_results data by celltype:
        cell_results_=cell_results[cell_results['celltypes']==celltype]
        
        # Prepare cell_results data by patient ids:
        patient_ids=list(set(cell_results_.patient_id))
        mean_degree_centrality=[]
        mean_average_clustering=[]
        mean_closeness_centrality=[]
        condition=[]
        df_by_patient_ids=pd.DataFrame()
        df_by_patient_ids['patient_ids']=patient_ids
        for patientId in patient_ids:
            patient_res=cell_results_[cell_results_['patient_id']==patientId]
            mean_degree_centrality.append(np.mean(list(patient_res['degree_centrality'])))
            mean_average_clustering.append(np.mean(list(patient_res['average_clustering'])))
            mean_closeness_centrality.append(np.mean(list(patient_res['closeness_centrality'])))
            condition.append(list(patient_res.condition)[0])
        df_by_patient_ids['mean_degree_centrality']=mean_degree_centrality
        df_by_patient_ids['mean_average_clustering']=mean_average_clustering
        df_by_patient_ids['mean_closeness_centrality']=mean_closeness_centrality
        df_by_patient_ids['condition']=condition
        df_by_patient_ids['celltype']=celltype
        cell_results_by_patient_id.append(df_by_patient_ids)
    cell_results_by_patient_id=pd.concat(cell_results_by_patient_id, axis=0)
    return cell_results_by_patient_id

def _abbreviate_strings_(list_of_strings, n_start, n_end):
    new_list_of_strings=[]
    for string in list_of_strings:
        if len(string)>(n_start+n_end):
            new_list_of_strings.append(string[:n_start]+string[len(string) - n_end:])
        else:
            new_list_of_strings.append(string)
    return new_list_of_strings
        
    
def _make_plot_name_from_column_(dataframe):
    patientID_condition_lookup=dict(zip(list(dataframe.patient_id), list(dataframe.condition)))
    n_start=3
    n_end=3
    patientID_abbreviatedID_lookup={}
    for patientID, disease_condition in patientID_condition_lookup.items():
        if len(patientID)>(n_start+n_end):
            # patient_id_abbreviated=patientID[:n_start]+patientID[len(patientID) - n_end:]+r" \bf{("+str(disease_condition)+")}"                       
            patient_id_abbreviated=patientID[:n_start]+patientID[len(patientID) - n_end:]+f" ({disease_condition})"                   
            patientID_abbreviatedID_lookup[patientID]=patient_id_abbreviated
        else:
            # patientID_abbreviatedID_lookup[patientID]=patientID+r"\bf{("+str(disease_condition)+")}"
            patientID_abbreviatedID_lookup[patientID]=patientID+f" ({disease_condition})"
    return patientID_abbreviatedID_lookup












