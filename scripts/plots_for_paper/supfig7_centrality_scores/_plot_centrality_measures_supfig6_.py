import seaborn as sns
sns.set_theme(style="whitegrid")
import itertools
from statannotations.Annotator import Annotator
from decimal import Decimal
import matplotlib.pyplot as plt

def _plot_centrality_measures_(p_values_cell_type_squidpy_centralityScores, data, celltype, conditions, scores, subplot_axis_id, legend=None, legend_loc=None, palette=None):
    disease_combinations=list(itertools.combinations(conditions, 2))
    var_name='centrality_measure'
    value_name='centrality_score'
    data=data.melt('condition', var_name=var_name, value_name=value_name)
    args = dict(x=var_name, y=value_name, data=data, hue="condition", hue_order=list(conditions), order=scores)
    pairs=[]
    pvals=[]
    for k in scores:
        for j in disease_combinations:
            pairs.append(((k, j[0]), (k, j[1])))
            if len(p_values_cell_type_squidpy_centralityScores[p_values_cell_type_squidpy_centralityScores['condition_1']==j[0]][p_values_cell_type_squidpy_centralityScores['condition_2']==j[1]][p_values_cell_type_squidpy_centralityScores['score']==k])==0:
                p_=str(p_values_cell_type_squidpy_centralityScores[p_values_cell_type_squidpy_centralityScores['condition_1']==j[1]][p_values_cell_type_squidpy_centralityScores['condition_2']==j[0]][p_values_cell_type_squidpy_centralityScores['score']==k]['p_value'][p_values_cell_type_squidpy_centralityScores['celltypes']==celltype].values[0])
            else:
                p_=str(p_values_cell_type_squidpy_centralityScores[p_values_cell_type_squidpy_centralityScores['condition_1']==j[0]][p_values_cell_type_squidpy_centralityScores['condition_2']==j[1]][p_values_cell_type_squidpy_centralityScores['score']==k]['p_value'][p_values_cell_type_squidpy_centralityScores['celltypes']==celltype].values[0])
            if float(p_)<0.05:
                pvals.append('%.2E' % Decimal(p_))
            else:
                pvals.append('$ns$')
    pairs_pvals_dict=dict(zip(pairs, pvals))
    fig, axes = plt.subplots(figsize=(20,10))
    sns.set(font_scale = 1.2)
    sns.set_style("white")
    ax = sns.violinplot(**args, cut=0, palette=palette)
    annot = Annotator(ax, pairs, **args)
    annot.set_custom_annotations(pvals)
    annot.annotate()
    plt.close()
    pairs_corrected=[]
    pvals_corrected=[]
    for j in range(len(annot.__dict__['annotations'])):
        pair_1=annot.__dict__['annotations'][j].__dict__['structs'][0]['group']
        pair_2=annot.__dict__['annotations'][j].__dict__['structs'][1]['group']   
        pairs_corrected.append((pair_1, pair_2))
        pvals_corrected.append(pairs_pvals_dict[(pair_1, pair_2)])
    
    plot_title=f'{celltype}'
    sns.set(font_scale = 1.2)
    sns.set_style("white")
    ax = sns.violinplot(ax=subplot_axis_id, **args, cut=0, palette=palette)
    annot = Annotator(ax, pairs, **args)
    annot.configure(text_format='simple', loc='inside', verbose=2) # fontsize=25
    annot.set_custom_annotations(pvals_corrected)
    annot.annotate()
    ax.set_xlabel('', fontsize=20, labelpad=25)
    ax.set_ylabel('Centrality scores', fontsize=25) # fontsize=25, labelpad=20
    ax.set_title(plot_title, fontsize=30) # fontsize=25, pad=20
    # sns.move_legend(
    #     ax, "lower right",
    #     bbox_to_anchor=(1.02, -0.04), ncol=3, title=None, frameon=False,
    # )
    yticks=[]
    for j in ax.get_yticks():
        yticks.append(round(j,1))
    ax.set_yticklabels(yticks, fontsize=20) # size = 20
    ax.set_xticklabels(['Degree centrality', 'Closeness centrality'], fontsize=20) # size = 20
    if legend:
        sns.move_legend(ax, loc=legend_loc,
                        title=None, fontsize=20,
                        bbox_to_anchor=(1.5, 0.92))
    else:
        ax.legend([],[], frameon=False)
    ax.grid(False)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
        
        