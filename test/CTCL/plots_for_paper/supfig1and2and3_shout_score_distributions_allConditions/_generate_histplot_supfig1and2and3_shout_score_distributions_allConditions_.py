import seaborn as sns
sns.set_theme(style="whitegrid")

def _generate_histplot_(score, data, subplot_axis_id, hue, legend_loc, title):
    # ax=sns.jointplot(ax=subplot_axis_id, data=data, x='local_entropy_2', y='egophily_2', hue='condition')
    ax=sns.histplot(ax=subplot_axis_id, data=data, x=score, hue=hue, kde=True)
    ax.set_title(title)
    sns.move_legend(ax, loc=legend_loc, title=None)