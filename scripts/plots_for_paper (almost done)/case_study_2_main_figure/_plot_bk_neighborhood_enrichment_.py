import seaborn as sns
sns.set_theme(style="whitegrid")
import itertools
from statannotations.Annotator import Annotator
from decimal import Decimal

def _plot_bk_neighborhood_enrichment_(p_values_cell_type_squidpy_nhoodEnrichment, squidpy_nhoodEnrichment_results, data, celltype, reference_celltype, conditions, subplot_axis_id, title_prefix=None, plot_title=None, title_loc=None, suptitle=None, ylabel='Neighborhood enrichment\n(z-scores)', xlabel=None, palette=None):
    disease_combinations=list(itertools.combinations(conditions, 2))
    pairs=disease_combinations.copy()
    pvals=[]
    for j in disease_combinations:
        if len(p_values_cell_type_squidpy_nhoodEnrichment[p_values_cell_type_squidpy_nhoodEnrichment['condition_1']==j[0]][p_values_cell_type_squidpy_nhoodEnrichment['condition_2']==j[1]])==0:
            p_df=p_values_cell_type_squidpy_nhoodEnrichment[p_values_cell_type_squidpy_nhoodEnrichment['condition_1']==j[1]][p_values_cell_type_squidpy_nhoodEnrichment['condition_2']==j[0]]
        else:
            p_df=p_values_cell_type_squidpy_nhoodEnrichment[p_values_cell_type_squidpy_nhoodEnrichment['condition_1']==j[0]][p_values_cell_type_squidpy_nhoodEnrichment['condition_2']==j[1]]
        
        if len(p_df[p_df['celltype_1']==celltype][p_df['celltype_2']==reference_celltype])==0:
            p_=str(p_df[p_df['celltype_1']==reference_celltype][p_df['celltype_2']==celltype]['p_value_adj'].values[0])
        else:
            p_=str(p_df[p_df['celltype_1']==celltype][p_df['celltype_2']==reference_celltype]['p_value_adj'].values[0])
        
        if float(p_)<0.05:
            pvals.append('%.2E' % Decimal(p_))
        else:
            pvals.append('$ns$')    
    if palette:
        args = dict(x="condition", y="sq_ne_zscore", data=data, order=list(conditions), palette=palette) # hue="condition", order=list(conditions_)
    else:
        args = dict(x="condition", y="sq_ne_zscore", data=data, order=list(conditions)) # hue="condition", order=list(conditions_)
    if title_prefix:
        plot_title_=r"$\bf{" + title_prefix + "}$"
        if suptitle:
            plot_title_+=suptitle
            plot_title=plot_title_+'\n'+plot_title
        else:
            plot_title=plot_title_+plot_title
    sns.set(font_scale = 1.2)
    sns.set_style("white")
    ax = sns.violinplot(ax=subplot_axis_id, **args, cut=0)
    ax.set_title(plot_title, loc=title_loc, fontsize=25) # pad=2
    ax.set_xlabel(xlabel, fontsize=25, ) # labelpad=5
    ax.set_ylabel(ylabel, fontsize=25, labelpad=5)
    annot = Annotator(ax, pairs, **args)
    annot.configure(text_format='simple', loc='inside', verbose=2, fontsize=20) # fontsize=25
    annot.set_custom_annotations(pvals)
    annot.annotate()
    yticks=[]
    for j in ax.get_yticks():
        yticks.append(round(j,1))
    ax.set_yticklabels(yticks, size = 20)
    ax.tick_params(axis='x', which='major', labelsize=20)
    ax.grid(False)
    
    
    