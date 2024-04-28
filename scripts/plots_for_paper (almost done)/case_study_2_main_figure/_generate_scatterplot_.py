import seaborn as sns
sns.set_theme(style="whitegrid")

def _generate_scatterplot_(x_axis_param, y_axis_param, data, subplot_axis_id, hue, legend_loc, title_prefix):
    # ax=sns.jointplot(ax=subplot_axis_id, data=data, x='local_entropy_2', y='egophily_2', hue='condition')
    # Note: jointplot (problem id: 1), lmlplot (p_id: 2), regplot(p_id: 2) 
    plot_title='Scatter plot of entropy and egophily scores\n$(radius=3)$'
    if title_prefix:
        plot_title=r"$\bf{" + title_prefix + "}$" + plot_title
    ax=sns.scatterplot(ax=subplot_axis_id, data=data, x=x_axis_param, y=y_axis_param, hue=hue) # s=10, alpha=0.25                                
    ax.set_title(plot_title)
    sns.move_legend(ax, loc=legend_loc, title=None)