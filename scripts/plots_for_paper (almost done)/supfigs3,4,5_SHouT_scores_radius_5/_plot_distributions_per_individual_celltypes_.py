import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
import itertools
from statannotations.Annotator import Annotator
from decimal import Decimal

def _plot_distributions_per_individual_celltypes_(p_values_cell_type, cell_results, celltype, heterogeneity_measure, conditions, radii, subplot_axis_id, legend_loc=None, title_prefix=None, plot_title=None, title_loc=None, xlabel=None, palette=None, legend=None):
    disease_combinations=list(itertools.combinations(conditions, 2))
    if heterogeneity_measure=='egophily':
        local_heterogeneity_measure=[f'{heterogeneity_measure}_{radius}' for radius in radii]
        var_name=f'{heterogeneity_measure}_measure'
        value_name=f'{heterogeneity_measure}_score'
        ylabel='Egophily'
    else:
        local_heterogeneity_measure=[f'local_{heterogeneity_measure}_{radius}' for radius in radii]
        var_name=f'local_{heterogeneity_measure}_measure'
        value_name=f'local_{heterogeneity_measure}_score'
        ylabel=f'Local {heterogeneity_measure}'
    columns=local_heterogeneity_measure+['condition']
    data_=cell_results[cell_results['cell_type']==celltype][columns]
    data=data_.melt('condition', var_name=var_name, value_name=value_name)
    if palette:
        args = dict(x=var_name, y=value_name, data=data, hue="condition", hue_order=list(conditions), order=local_heterogeneity_measure, palette=palette)
    else:
        args = dict(x=var_name, y=value_name, data=data, hue="condition", hue_order=list(conditions), order=local_heterogeneity_measure)
    pairs=[]
    pvals=[]
    for j in local_heterogeneity_measure:
        for k in disease_combinations:
            pairs.append(((j, k[0]), (j, k[1])))
            if len(p_values_cell_type[p_values_cell_type['cell_type']==celltype][p_values_cell_type['condition_1']==k[0]][p_values_cell_type['condition_2']==k[1]][p_values_cell_type['score']==j]['p_value_adj'])==0:
                p_=str(p_values_cell_type[p_values_cell_type['cell_type']==celltype][p_values_cell_type['condition_1']==k[1]][p_values_cell_type['condition_2']==k[0]][p_values_cell_type['score']==j]['p_value_adj'].values[0])
            else:
                p_=str(p_values_cell_type[p_values_cell_type['cell_type']==celltype][p_values_cell_type['condition_1']==k[0]][p_values_cell_type['condition_2']==k[1]][p_values_cell_type['score']==j]['p_value_adj'].values[0])
            if float(p_)<0.05:
                pvals.append('%.2E' % Decimal(p_))
            else:
                pvals.append('$ns$')
    pairs_pvals_dict=dict(zip(pairs, pvals))
    fig, axes = plt.subplots(figsize=(20,10))
    sns.set(font_scale = 1.2)
    sns.set_style("white")
    ax = sns.violinplot(**args, cut=0)
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
    sns.set(font_scale = 1.2)
    sns.set_style("white")
    ax = sns.violinplot(ax=subplot_axis_id, **args, cut=0)
    if title_prefix:
        if plot_title:
            plot_title=r"$\bf{" + title_prefix + "}$" + plot_title
        else:
            plot_title=r"$\bf{" + title_prefix + "}$"
    ax.set_title(plot_title, loc=title_loc, fontsize=25) # fontsize=25, pad=20
    ax.set_xlabel(xlabel, fontsize=25) # fontsize=25, labelpad=20
    ax.set_ylabel(ylabel, fontsize=25, labelpad=5) # fontsize=25, labelpad=20
    ax.set(xticklabels=[]) # size=20
    annot = Annotator(ax, pairs, **args)
    annot.configure(text_format='simple', loc='inside', verbose=2, fontsize=20) # fontsize=25
    annot.set_custom_annotations(pvals_corrected)
    annot.annotate()
    # ax.legend(title='Condition', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    # ax.legend(title=None, loc=legend_loc, borderaxespad=0, fontsize=20)
    ax.get_legend().set_visible(False)
    # ax.legend.remove()
    yticks=[]
    for j in ax.get_yticks():
        yticks.append(round(j,1))
    ax.set_yticklabels(yticks, size = 20) # size = 20
    ax.tick_params(axis='x', which='major', labelsize=20)
    ax.grid(False)
    if legend:
        sns.move_legend(ax, loc=legend_loc,
                        title=None, fontsize=20,
                        bbox_to_anchor=(1.5, 0.85))
    else:
        ax.legend([],[], frameon=False)
    
    
    
    