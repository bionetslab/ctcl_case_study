import pandas as pd
import scipy.stats as stats
import os
import itertools as itt
from IPython.display import set_matplotlib_formats
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import itertools
from statannotations.Annotator import Annotator
from decimal import Decimal
import numpy as np
from decimal import Decimal

def plot_heatmap(data, genes, figsize=(9,6), outfile=None, title=None, subplot_axis_id=None):
    if subplot_axis_id:
        ax=sns.heatmap(data=data[genes],ax=subplot_axis_id)
    else:
        fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(data=data[genes],ax=ax)
    ax.set(ylabel='')
    if not title is None:
        ax.set_title(title)
    # fig.tight_layout()
    if not outfile is None:
        fig.savefig(outfile)

if __name__ == '__main__':
    p_values_cell_type=pd.read_csv(os.path.join('results', 'p_values_cell_type.csv'))
    p_values_cell_type_shuffled_labels=pd.read_csv(os.path.join('pvalue_statistics_shuffled_labels.csv'))
    radii=[1,2,3,4,5]
    het_measures=['entropy', 'egophily', 'homophily']
    # layout=[
    #     [het_measures[0], het_measures[0], het_measures[1], het_measures[1]],
    #     [".", het_measures[2], het_measures[2], "."]
    # ]
    # layout=[]
    # for hmeasure in het_measures:
    #     layout.append([hmeasure])
    celltypes=list(np.unique(p_values_cell_type_shuffled_labels['cell_type']))
    conditions_=list(set(list(np.unique(p_values_cell_type_shuffled_labels.condition_1))).union(set(list(np.unique(p_values_cell_type_shuffled_labels.condition_2)))))
    conditions_abbreviations_dict=dict(zip(conditions_, ['AD', 'PSO', 'CTCL']))
    disease_combinations=list(itertools.combinations(conditions_, 2))
    layout=[
        [str(disease_combinations[0]), str(disease_combinations[0]), str(disease_combinations[1]), str(disease_combinations[1])],
        [".", str(disease_combinations[2]), str(disease_combinations[2]), "."]
    ]
    for het_measure in het_measures:
        if het_measure=='egophily':
            local_het_measure=f'{het_measure}'
        else:
            local_het_measure=f'local_{het_measure}'
        fig, axes = plt.subplot_mosaic(layout, layout='constrained', figsize=(20,10))
        local_het_measures=[]
        pvalue_data_all_radii=[]
        pvalue_actual_data_all_radii=[]
        for radius in radii:
            het_measure_name=f'{local_het_measure}_{radius}'
            local_het_measures.append(het_measure_name)
            pvalue_data_per_radius=p_values_cell_type_shuffled_labels[p_values_cell_type_shuffled_labels['score']==het_measure_name]
            pvalue_actual_data_per_radius=p_values_cell_type[p_values_cell_type['score']==het_measure_name]
            pvalue_data_all_radii.append(pvalue_data_per_radius)
            pvalue_actual_data_all_radii.append(pvalue_actual_data_per_radius)
        pvalue_data_all_radii=pd.concat(pvalue_data_all_radii, axis=0)
        pvalue_actual_data_all_radii=pd.concat(pvalue_actual_data_all_radii, axis=0)
        for disease_combination in disease_combinations:
            # ----- Make subplot: -----
            subplot_axis_id=axes[str(disease_combination)]
            legend_loc='lower right'
            title_prefix='A. '
            # -------------------------
            if len(pvalue_data_all_radii[pvalue_data_all_radii['condition_1']==disease_combination[0]][pvalue_data_all_radii['condition_2']==disease_combination[1]])==0:
                data_for_heatmap_subplot=pvalue_data_all_radii[pvalue_data_all_radii['condition_1']==disease_combination[1]][pvalue_data_all_radii['condition_2']==disease_combination[0]]
            else:
                data_for_heatmap_subplot=pvalue_data_all_radii[pvalue_data_all_radii['condition_1']==disease_combination[0]][pvalue_data_all_radii['condition_2']==disease_combination[1]]
            actual_data_for_heatmap_subplot=pvalue_actual_data_all_radii[pvalue_actual_data_all_radii['condition_1']==disease_combination[1]][pvalue_actual_data_all_radii['condition_2']==disease_combination[0]]
            if len(actual_data_for_heatmap_subplot)==0: # The (condition_1, condition_2) column pair values are just filled (because the snytzetic data contains pairs in alphabetical order, i.e., condition_1 value always smaller than condition_2 value alphabetically.)
                actual_data_for_heatmap_subplot=pvalue_actual_data_all_radii[pvalue_actual_data_all_radii['condition_1']==disease_combination[0]][pvalue_actual_data_all_radii['condition_2']==disease_combination[1]]
            # Generate local p-value statistic matrices for each subplot:
            p_value_mean=[]
            actual_data_p_value=[]
            columns=list(set(data_for_heatmap_subplot['score']))
            for celltype in celltypes:
                data_for_heatmap_subplot_filtered_by_celltype=data_for_heatmap_subplot[data_for_heatmap_subplot['cell_type']==celltype]
                actual_data_for_heatmap_subplot_filtered_by_celltype=actual_data_for_heatmap_subplot[actual_data_for_heatmap_subplot['cell_type']==celltype]
                df_p_value_mean=pd.DataFrame()
                df_p_value_mean['cell_type']=[celltype]
                df_actual_data_p_value=pd.DataFrame()
                df_actual_data_p_value['cell_type']=[celltype]
                mean_scores_dict={'cell_type': celltype}
                actual_data_scores_dict={'cell_type': celltype}
                for score in columns:
                    mean_scores_dict[score]=[data_for_heatmap_subplot_filtered_by_celltype[data_for_heatmap_subplot_filtered_by_celltype['score']==score]['p_value_mean'].values[0]]
                    actual_data_scores_dict[score]=[actual_data_for_heatmap_subplot_filtered_by_celltype[actual_data_for_heatmap_subplot_filtered_by_celltype['score']==score]['p_value'].values[0]]
                mean_scores_df=pd.DataFrame.from_dict(mean_scores_dict)
                actual_data_scores_df=pd.DataFrame.from_dict(actual_data_scores_dict)
                p_value_mean.append(mean_scores_df)
                actual_data_p_value.append(actual_data_scores_df)
            p_value_mean=pd.concat(p_value_mean, axis=0).set_index('cell_type', drop=True)
            actual_data_p_value=pd.concat(actual_data_p_value, axis=0).set_index('cell_type', drop=True)
            robustness_score=p_value_mean-actual_data_p_value # Positive values are preferred. Higher the value, more robust the SHouT algorithm for generating local heterogeneity scores.
            plot_heatmap(robustness_score, 
                         local_het_measures,
                         title=str(disease_combination),
                         subplot_axis_id=subplot_axis_id)
        # ========== Generate, save and show final plot (fig1): ==========
        plt.suptitle(f'{local_het_measure}')
        plt.savefig(f'robustnessScores_shuffledLabels_pvalueMean_{local_het_measure}.pdf', format='pdf', bbox_inches='tight') # Figure name options
        plt.show()
    
    