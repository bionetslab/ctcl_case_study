import seaborn as sns
sns.set_theme(style="whitegrid")
import itertools
from statannotations.Annotator import Annotator
from decimal import Decimal

def _plot_bk_neighborhood_enrichment_(p_values_cell_type_squidpy_nhoodEnrichment, squidpy_nhoodEnrichment_results, data, celltype, conditions, radii, subplot_axis_id, title_prefix, suptitle):
    disease_combinations=list(itertools.combinations(conditions, 2))
    pairs=disease_combinations.copy()
    pvals=[]
    for j in disease_combinations:
        if len(p_values_cell_type_squidpy_nhoodEnrichment[p_values_cell_type_squidpy_nhoodEnrichment['condition_1']==j[0]][p_values_cell_type_squidpy_nhoodEnrichment['condition_2']==j[1]])==0:
            p_df=p_values_cell_type_squidpy_nhoodEnrichment[p_values_cell_type_squidpy_nhoodEnrichment['condition_1']==j[1]][p_values_cell_type_squidpy_nhoodEnrichment['condition_2']==j[0]]
        else:
            p_df=p_values_cell_type_squidpy_nhoodEnrichment[p_values_cell_type_squidpy_nhoodEnrichment['condition_1']==j[0]][p_values_cell_type_squidpy_nhoodEnrichment['condition_2']==j[1]]
        
        if len(p_df[p_df['celltype_1']==celltype][p_df['celltype_2']=='Basal keratinocytes'])==0:
            p_=str(p_df[p_df['celltype_1']=='Basal keratinocytes'][p_df['celltype_2']==celltype]['p_value'].values[0])
        else:
            p_=str(p_df[p_df['celltype_1']==celltype][p_df['celltype_2']=='Basal keratinocytes']['p_value'].values[0])
        
        if float(p_)<0.05:
            pvals.append('%.2E' % Decimal(p_))
        else:
            pvals.append('$ns$')    
    args = dict(x="condition", y="sq_ne_zscore", data=data, order=list(conditions)) # hue="condition", order=list(conditions_)
    plot_title=str(f'{celltype}')
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
    ax.set_title(plot_title, ) # fontsize=25, pad=20
    ax.set_xlabel('radius', ) # fontsize=25, labelpad=20
    ax.set_ylabel('nhood enrichment', ) # fontsize=25, labelpad=20
    annot = Annotator(ax, pairs, **args)
    annot.configure(text_format='simple', loc='inside', verbose=2, ) # fontsize=25
    annot.set_custom_annotations(pvals)
    annot.annotate()
    xticks=[]
    for j in ax.get_yticks():
        xticks.append(round(j,1))
    ax.set_yticklabels(xticks, ) # size = 20
    
    
    
    
    