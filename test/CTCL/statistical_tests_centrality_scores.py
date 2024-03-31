import pandas as pd
import scipy.stats as stats
import os
import itertools as itt

if __name__ == '__main__':
    p_values_global = {
        'condition_1': [],
        'condition_2': [],
        'score': [],
        'statistic': [],
        'p_value': [],
        'p_value_adj': []
    }
    p_values_cell_type = {
        'condition_1': [],
        'condition_2': [],
        'score': [],
        'celltypes': [],
        'statistic': [],
        'p_value': [],
        'p_value_adj': []
    }
    scores = ['degree_centrality', 'average_clustering', 'closeness_centrality']
    squidpy_nhoodEnrichment_results = pd.read_csv(os.path.join('results', 'squidpy_centralityScores_results.csv'))
    conditions = list(set(squidpy_nhoodEnrichment_results.condition))
    cell_types = list(set(squidpy_nhoodEnrichment_results.celltypes))
    filters = {condition: (squidpy_nhoodEnrichment_results.condition == condition) for condition in conditions}
    for condition_1, condition_2 in itt.combinations(conditions, 2):
        for score in scores:
            column=scores+['celltypes']
            scores_1 = squidpy_nhoodEnrichment_results[column][filters[condition_1]]
            scores_2 = squidpy_nhoodEnrichment_results[column][filters[condition_2]]
            statistic, p_value = stats.mannwhitneyu(x=scores_1[score], y=scores_2[score])
            p_values_global['condition_1'].append(condition_1)
            p_values_global['condition_2'].append(condition_2)
            p_values_global['score'].append(score)
            p_values_global['statistic'].append(statistic)
            p_values_global['p_value'].append(p_value)
            p_values_global['p_value_adj'].append(p_value * 3)
            for celltype in cell_types:
                results_condition_1_celltype = scores_1[scores_1.celltypes == celltype]
                results_condition_2_celltype = scores_2[scores_2.celltypes == celltype]
                scores_1_celltype = results_condition_1_celltype[score]
                scores_2_celltype = results_condition_2_celltype[score]
                statistic, p_value = stats.mannwhitneyu(x=scores_1_celltype, y=scores_2_celltype)
                p_values_cell_type['condition_1'].append(condition_1)
                p_values_cell_type['condition_2'].append(condition_2)
                p_values_cell_type['score'].append(score)
                p_values_cell_type['celltypes'].append(celltype)
                p_values_cell_type['statistic'].append(statistic)
                p_values_cell_type['p_value'].append(p_value)
                p_values_cell_type['p_value_adj'].append(p_value * 3 * len(cell_types) * 2)
        
    p_values_global = pd.DataFrame(data=p_values_global)
    p_values_cell_type = pd.DataFrame(data=p_values_cell_type)
    p_values_global.to_csv(os.path.join('results', 'p_values_global_squidpy_centralityScores.csv'))
    p_values_cell_type.to_csv(os.path.join('results', 'p_values_cell_type_squidpy_centralityScores.csv'))
