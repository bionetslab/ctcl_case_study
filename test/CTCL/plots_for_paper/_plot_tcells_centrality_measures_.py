import seaborn as sns
sns.set_theme(style="whitegrid")
import itertools
from statannotations.Annotator import Annotator
from decimal import Decimal

def _plot_tcells_centrality_measures_(p_values_cell_type_squidpy_centralityScores, data, celltype, conditions, radii, scores, subplot_axis_id, title_prefix):
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
    plot_title=f'Centrality scores ({celltype})'
    if title_prefix:
        plot_title=r"$\bf{" + title_prefix + "}$" + plot_title
    sns.set(font_scale = 1.2)
    sns.set_style("white")
    ax = sns.violinplot(ax=subplot_axis_id, **args, cut=0)
    annot = Annotator(ax, pairs, **args)
    annot.configure(text_format='simple', loc='inside', verbose=2, ) # fontsize=25
    annot.set_custom_annotations(pvals)
    annot.annotate()
    ax.set_xlabel('', fontsize=25, labelpad=20)
    ax.set_ylabel('Centrality scores', ) # fontsize=25, labelpad=20
    ax.set_title(plot_title, ) # fontsize=25, pad=20
    sns.move_legend(
        ax, "lower right",
        bbox_to_anchor=(1.03, -0.08), ncol=3, title=None, frameon=False,
    )
    xticks=[]
    for j in ax.get_yticks():
        xticks.append(round(j,1))
    ax.set_yticklabels(xticks, ) # size = 20
        
        
        