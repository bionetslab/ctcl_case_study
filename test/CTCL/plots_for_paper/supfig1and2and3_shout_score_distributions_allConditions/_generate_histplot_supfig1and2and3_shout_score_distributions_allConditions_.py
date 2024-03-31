import seaborn as sns
sns.set_theme(style="whitegrid")

def _generate_histplot_(score, data, subplot_axis_id, hue, legend_loc, title, legend=None):
    ax=sns.histplot(ax=subplot_axis_id, data=data, x=score, hue=hue, kde=True)
    if legend:
        sns.move_legend(ax, loc=legend_loc,
                        title=None, fontsize=20,
                        bbox_to_anchor=(1.5, 0.5))
    else:
        ax.legend([],[], frameon=False)
    ax.set_title(title, loc='left', fontsize=25) # fontsize=25, pad=20
    ax.set_xlabel(score, fontsize=20) # fontsize=25, labelpad=20
    ax.set_ylabel('Count', fontsize=20, labelpad=5) # fontsize=25, labelpad=20
    yticks=[]
    for j in ax.get_yticks():
        yticks.append(int(j))
    ax.set_yticklabels(yticks, size = 20) # size = 20
    ax.tick_params(axis='x', which='major', labelsize=20)
    ax.grid(False)