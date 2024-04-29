import pandas as pd
import scipy.stats as stats
import os
import itertools as itt

if __name__ == '__main__':
    N=15
    p_values_global_subsampled_patientwise=[]
    p_values_cell_type_subsampled_patientwise=[]
    for i in range(100):
        print(i)
        p_values_global = {
            'condition_1': [],
            'condition_2': [],
            'score': [],
            'statistic': [],
            'p_value': [],
            'p_value_adj': []
        }
        global_scores = ['global_homophily', 'global_entropy']
        sample_results = pd.read_csv(os.path.join('../../results', 'sample_results.csv'))
        conditions = list(set(sample_results.condition))
        filters = {condition: (sample_results.condition == condition) for condition in conditions}
        patient_ids=[]
        for condition_1, condition_2 in itt.combinations(conditions, 2):
            sample_results_condition_1=sample_results[sample_results['condition']==condition_1].copy().sample(n=N)
            sample_results_condition_2=sample_results[sample_results['condition']==condition_2].copy().sample(n=N)
            sample_results_subsampled_patientwise=pd.concat([sample_results_condition_1, sample_results_condition_2], axis=0)
            patient_ids_=list(sample_results_subsampled_patientwise.patient_id)
            patient_ids.append(patient_ids_)
            for score in global_scores:
                scores_1 = sample_results_subsampled_patientwise[score][filters[condition_1]]
                scores_2 = sample_results_subsampled_patientwise[score][filters[condition_2]]
                statistic, p_value = stats.mannwhitneyu(x=scores_1, y=scores_2)
                p_values_global['condition_1'].append(condition_1)
                p_values_global['condition_2'].append(condition_2)
                p_values_global['score'].append(score)
                p_values_global['statistic'].append(statistic)
                p_values_global['p_value'].append(p_value)
                p_values_global['p_value_adj'].append(p_value * 3)
        p_values_cell_type = {
            'condition_1': [],
            'condition_2': [],
            'score': [],
            'cell_type': [],
            'statistic': [],
            'p_value': [],
            'p_value_adj': []
        }
        score_types = ['local_homophily', 'local_entropy', 'egophily']
        local_scores = [f'{score_type}_{radius}' for score_type in score_types for radius in [1, 2, 3, 4, 5]]
        cell_results = pd.read_csv(os.path.join('../../results', 'cell_results.csv'))
        cnt=-1
        for condition_1, condition_2 in itt.combinations(conditions, 2):
            cnt+=1
            cell_results_=[]
            for patientId in patient_ids[cnt]:
                cell_results_.append(cell_results[cell_results['patient_id']==patientId])
            cell_results_=pd.concat(cell_results_, axis=0)
            filters = {condition: (cell_results_.condition == condition) for condition in conditions}
            cell_types = list(set(cell_results_.cell_type))
            for score in local_scores:
                results_condition_1 = cell_results_[filters[condition_1]]
                results_condition_2 = cell_results_[filters[condition_2]]
                scores_1 = results_condition_1[score]
                scores_2 = results_condition_2[score]
                statistic, p_value = stats.mannwhitneyu(x=scores_1, y=scores_2)
                p_values_global['condition_1'].append(condition_1)
                p_values_global['condition_2'].append(condition_2)
                p_values_global['score'].append(score)
                p_values_global['statistic'].append(statistic)
                p_values_global['p_value'].append(p_value)
                p_values_global['p_value_adj'].append(p_value * 3)
                for cell_type in cell_types:
                    results_condition_1_celltype = results_condition_1[results_condition_1.cell_type == cell_type]
                    results_condition_2_celltype = results_condition_2[results_condition_2.cell_type == cell_type]
                    scores_1 = results_condition_1_celltype[score]
                    scores_2 = results_condition_2_celltype[score]
                    statistic, p_value = stats.mannwhitneyu(x=scores_1, y=scores_2)
                    p_values_cell_type['condition_1'].append(condition_1)
                    p_values_cell_type['condition_2'].append(condition_2)
                    p_values_cell_type['score'].append(score)
                    p_values_cell_type['cell_type'].append(cell_type)
                    p_values_cell_type['statistic'].append(statistic)
                    p_values_cell_type['p_value'].append(p_value)
                    p_values_cell_type['p_value_adj'].append(p_value * 3 * len(cell_types) * 5)
        p_values_global = pd.DataFrame(data=p_values_global)
        p_values_cell_type = pd.DataFrame(data=p_values_cell_type)
        p_values_global_subsampled_patientwise.append(p_values_global)
        p_values_cell_type_subsampled_patientwise.append(p_values_cell_type)
    p_values_global_subsampled_patientwise=pd.concat(p_values_global_subsampled_patientwise, axis=0)
    p_values_cell_type_subsampled_patientwise=pd.concat(p_values_cell_type_subsampled_patientwise, axis=0)
    p_values_global_subsampled_patientwise.to_csv(os.path.join('../../results', 'p_values_global_subsampled_patientwise.csv'))
    p_values_cell_type_subsampled_patientwise.to_csv(os.path.join('../../results', 'p_values_cell_type_subsampled_patientwise.csv'))
