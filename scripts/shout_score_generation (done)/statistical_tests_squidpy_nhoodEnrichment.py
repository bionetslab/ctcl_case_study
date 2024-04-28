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
        'celltype_1': [],
        'celltype_2': [],
        'statistic': [],
        'p_value': [],
        'p_value_adj': []
    }
    score = 'sq_ne_zscore'
    squidpy_nhoodEnrichment_results = pd.read_csv(os.path.join('../../results', 'squidpy_nhoodEnrichment_results.csv'))
    conditions = list(set(squidpy_nhoodEnrichment_results.condition))
    cell_types = list(set(squidpy_nhoodEnrichment_results.celltype_1).union(set(squidpy_nhoodEnrichment_results.celltype_2)))
    filters = {condition: (squidpy_nhoodEnrichment_results.condition == condition) for condition in conditions}
    for condition_1, condition_2 in itt.combinations(conditions, 2):
        scores_1 = squidpy_nhoodEnrichment_results[[score, 'celltype_1', 'celltype_2']][filters[condition_1]]
        scores_2 = squidpy_nhoodEnrichment_results[[score, 'celltype_1', 'celltype_2']][filters[condition_2]]
        statistic, p_value = stats.mannwhitneyu(x=scores_1[score], y=scores_2[score])
        p_values_global['condition_1'].append(condition_1)
        p_values_global['condition_2'].append(condition_2)
        p_values_global['score'].append(score)
        p_values_global['statistic'].append(statistic)
        p_values_global['p_value'].append(p_value)
        p_values_global['p_value_adj'].append(p_value * 3)
        for celltype_1, celltype_2 in itt.combinations(cell_types, 2):
            celltype_1_, celltype_2_ = sorted([celltype_1, celltype_2])
            results_condition_1_celltype_1 = scores_1[scores_1.celltype_1 == celltype_1_][scores_1.celltype_2 == celltype_2_]
            results_condition_1_celltype_2 = scores_1[scores_1.celltype_2 == celltype_1_][scores_1.celltype_1 == celltype_2_]
            results_condition_1_celltype=pd.concat([results_condition_1_celltype_1, results_condition_1_celltype_2], axis=0)
            results_condition_2_celltype_1 = scores_2[scores_2.celltype_1 == celltype_1_][scores_2.celltype_2 == celltype_2_]
            results_condition_2_celltype_2 = scores_2[scores_2.celltype_2 == celltype_1_][scores_2.celltype_1 == celltype_2_]
            results_condition_2_celltype=pd.concat([results_condition_2_celltype_1, results_condition_2_celltype_2], axis=0)
            results_condition_2_celltype = scores_2[scores_2.celltype_2 == celltype_2_]
            scores_1_celltype = results_condition_1_celltype[score]
            scores_2_celltype = results_condition_2_celltype[score]
            statistic, p_value = stats.mannwhitneyu(x=scores_1_celltype, y=scores_2_celltype)
            p_values_cell_type['condition_1'].append(condition_1)
            p_values_cell_type['condition_2'].append(condition_2)
            p_values_cell_type['score'].append(score)
            p_values_cell_type['celltype_1'].append(celltype_1_)
            p_values_cell_type['celltype_2'].append(celltype_2_)
            p_values_cell_type['statistic'].append(statistic)
            p_values_cell_type['p_value'].append(p_value)
            p_values_cell_type['p_value_adj'].append(p_value * 3 * len(cell_types) * len(cell_types))
        
    p_values_global = pd.DataFrame(data=p_values_global)
    p_values_cell_type = pd.DataFrame(data=p_values_cell_type)
    p_values_global.to_csv(os.path.join('../../results', 'p_values_global_squidpy_nhoodEnrichment.csv'))
    p_values_cell_type.to_csv(os.path.join('../../results', 'p_values_cell_type_squidpy_nhoodEnrichment.csv'))
