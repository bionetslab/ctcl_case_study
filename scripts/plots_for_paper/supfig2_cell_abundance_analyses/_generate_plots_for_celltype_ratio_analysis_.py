import itertools as itt
from decimal import Decimal
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator

def _generate_violinplot_(data, categories, x_col_name, y_col_name, p_value, subplot_axis_id, legend_loc, title_prefix=None, plot_title='', pval_col1='condition_1', pval_col2='condition_2', xlabel=None, ylabel='y', order=None, palette=None):
    args = dict(x=x_col_name, y=y_col_name, data=data, order=categories)
    if order:
        args = dict(x=x_col_name, y=y_col_name, data=data, order=order)
    category_combinations=list(itt.combinations(categories, 2))
    category_tuple_orders_in_pvalue_data=list(zip(list(p_value.condition_1), list(p_value.condition_2)))
    pairs=[]
    pvals=[]
    
    for k in category_combinations:
        if k in category_tuple_orders_in_pvalue_data:
            pairs.append((k[0], k[1]))
            p_=str(p_value[p_value['condition_1']==k[0]][p_value['condition_2']==k[1]]['p'].values[0])
        else:
            pairs.append((k[0], k[1]))
            p_=str(p_value[p_value['condition_1']==k[1]][p_value['condition_2']==k[0]]['p'].values[0])
        if float(p_)<0.05:
            pvals.append('%.2E' % Decimal(p_))
        else:
            pvals.append('$ns$')
    pairs_pvals_dict=dict(zip(pairs, pvals))
    fig, axes = plt.subplots()
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
        try:
            pvals_corrected.append(pairs_pvals_dict[(pair_1, pair_2)])
        except:
            pvals_corrected.append(pairs_pvals_dict[(pair_2, pair_1)])
    sns.set(font_scale = 1.2)
    sns.set_style("white")
    ax = sns.violinplot(ax=subplot_axis_id, **args, cut=0, palette=palette)
    if title_prefix:
        plot_title=r"$\bf{" + title_prefix + "}$ " + plot_title
    ax.set_title(plot_title, fontsize=35) # fontsize=25, pad=20
    ax.set_xlabel(xlabel, fontsize=30) # fontsize=25, labelpad=20
    ax.set_ylabel(ylabel, fontsize=30) # fontsize=25, labelpad=20
    # ax.set_xticklabels(categories, ) # size=20
    annot = Annotator(ax, pairs, **args)
    annot.configure(text_format='simple', loc='inside', verbose=2, fontsize=30) # fontsize=25
    annot.set_custom_annotations(pvals_corrected)
    annot.annotate()
    # ax.legend(title='Condition', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    # ax.legend(title=None, loc=legend_loc, borderaxespad=0)
    yticks=[]
    for j in ax.get_yticks():
        yticks.append(round(j,1))
    ax.set_yticklabels(yticks, size = 30)
    ax.tick_params(axis='x', which='major', labelsize=30)
    ax.grid(False)


def _generate_histplot_(score, data, subplot_axis_id, hue, legend_loc, title_prefix=None, plot_title='', xlabel='Score', ylabel='Count', palette=None):
    ax=sns.histplot(ax=subplot_axis_id, data=data, x=score, hue=hue, kde=True, palette=palette)
    if title_prefix:
        plot_title=r"$\bf{" + title_prefix + "}$" + plot_title
    ax.set_title(plot_title) # fontsize=25, pad=20
    sns.move_legend(ax, loc=legend_loc, title=None)
    ax.grid(False)
    
    
    
    