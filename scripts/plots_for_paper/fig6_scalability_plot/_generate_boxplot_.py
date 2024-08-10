import seaborn as sns
sns.set_theme(style="whitegrid")

def _generate_boxplot_(data, var_name, val_name, subplot_axis_id, radii, legend_loc=None, title_prefix=None, plot_title=None, ylabel=None, palette=None, title_loc=None):
    args = dict(x=var_name, y=val_name, data=data, hue=var_name, hue_order=radii, palette=palette)
    # xticklabels=[f'$r={radius}$' for radius in radii]
    xticklabels=radii
    ax = sns.boxplot(**args, ax=subplot_axis_id, dodge=False)
    if title_prefix:
        if plot_title:
            plot_title=r"$\bf{" + title_prefix + "}$" + plot_title
        else:
            plot_title=r"$\bf{" + title_prefix + "}$"
    ax.set_title(plot_title, loc=title_loc, fontsize=25) # fontsize=25, pad=20
    # ax.set_xlabel(xlabel, fontsize=25) # fontsize=25, labelpad=20
    ax.set_xlabel('Radius ($r$)', fontsize=25, labelpad=5) # fontsize=25, labelpad=20
    ax.set_ylabel('Runtime, all samples (sec)', fontsize=25, labelpad=5) # fontsize=25, labelpad=20
    ax.set_xticklabels(xticklabels, size=20) # size=20
    
    ax.get_legend().set_visible(False)
    yticks=[]
    for j in ax.get_yticks():
        yticks.append(round(j,1))
    ax.set_yticklabels(yticks, size = 20) # size = 20
    ax.tick_params(axis='x', which='major', labelsize=20)
    ax.yaxis.set_tick_params(which='both', labelleft=True)
    ax.grid(False)
    
    
    
    