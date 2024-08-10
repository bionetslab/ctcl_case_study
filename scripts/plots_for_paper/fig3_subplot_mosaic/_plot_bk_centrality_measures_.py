import seaborn as sns
sns.set_theme(style="whitegrid")
import itertools
from statannotations.Annotator import Annotator
from decimal import Decimal

def _plot_bk_centrality_measures_(p_values_cell_type_squidpy_centralityScores, data, celltype, conditions, scores, subplot_axis_id, legend_loc=None, title_prefix=None, plot_title=None, title_loc=None, xlabel=None, ylabel=None, palette=None, ncol=None, bbox_to_anchor=None):
    disease_combinations=list(itertools.combinations(conditions, 2))
    var_name='centrality_measure'
    value_name='centrality_score'
    data=data.melt('condition', var_name=var_name, value_name=value_name)
    if palette:
        args = dict(x=var_name, y=value_name, data=data, hue="condition", hue_order=list(conditions), order=scores, palette=palette)
    else:
        args = dict(x=var_name, y=value_name, data=data, hue="condition", hue_order=list(conditions), order=scores)    
    pairs=[]
    pvals=[]
    for k in scores:
        for j in disease_combinations:
            pairs.append(((k, j[0]), (k, j[1])))
            if len(p_values_cell_type_squidpy_centralityScores[p_values_cell_type_squidpy_centralityScores['condition_1']==j[0]][p_values_cell_type_squidpy_centralityScores['condition_2']==j[1]][p_values_cell_type_squidpy_centralityScores['score']==k])==0:
                p_=str(p_values_cell_type_squidpy_centralityScores[p_values_cell_type_squidpy_centralityScores['condition_1']==j[1]][p_values_cell_type_squidpy_centralityScores['condition_2']==j[0]][p_values_cell_type_squidpy_centralityScores['score']==k]['p_value_adj'][p_values_cell_type_squidpy_centralityScores['celltypes']==celltype].values[0])
            else:
                p_=str(p_values_cell_type_squidpy_centralityScores[p_values_cell_type_squidpy_centralityScores['condition_1']==j[0]][p_values_cell_type_squidpy_centralityScores['condition_2']==j[1]][p_values_cell_type_squidpy_centralityScores['score']==k]['p_value_adj'][p_values_cell_type_squidpy_centralityScores['celltypes']==celltype].values[0])
            if float(p_)<0.05:
                pvals.append('%.2E' % Decimal(p_))
            else:
                pvals.append('$ns$')
    if title_prefix:
        if plot_title:
            plot_title=r"$\bf{" + title_prefix + "}$" + plot_title
        else:
            plot_title=r"$\bf{" + title_prefix + "}$"
    sns.set(font_scale = 1.2)
    sns.set_style("white")
    ax = sns.violinplot(ax=subplot_axis_id, **args, cut=0)
    annot = Annotator(ax, pairs, **args)
    annot.configure(text_format='simple', loc='inside', verbose=2, fontsize=20) # fontsize=25
    annot.set_custom_annotations(pvals)
    annot.annotate()
    ax.set_title(plot_title, loc=title_loc, fontsize=25) # fontsize=25, pad=20
    ax.set_xlabel(xlabel, fontsize=25) # fontsize=25, labelpad=20
    ax.set_ylabel(ylabel, fontsize=25, labelpad=5) # fontsize=25, labelpad=20
    if ncol:    
        sns.move_legend(
            ax, 
            loc=legend_loc,
            bbox_to_anchor=bbox_to_anchor,
            ncol=ncol,
            title=None,
            frameon=True,
            borderaxespad=0,
            fontsize=20
        )
    else:    
        sns.move_legend(
            ax, 
            loc=legend_loc,
            bbox_to_anchor=bbox_to_anchor,
            title=None,
            frameon=True,
            borderaxespad=0,
            fontsize=20
        )
    yticks=[]
    for j in ax.get_yticks():
        yticks.append(round(j,1))
    ax.set_yticklabels(yticks, size = 20) # size = 20
    # ax.tick_params(axis='x', which='major', labelsize=20) # or:
    ax.set_xticklabels(['Degree centrality', 'Closeness centrality'], size = 20) # size = 20
    ax.get_legend().remove()
    ax.grid(False)
        
        