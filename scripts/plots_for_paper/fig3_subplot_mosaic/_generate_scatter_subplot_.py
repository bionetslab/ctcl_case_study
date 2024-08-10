import seaborn as sns

def _generate_scatter_subplot_(background_image, scatter_plot_data, x, y, celltype_column_name, axes, palette='Set3', ncol=3, title_prefix=None, plot_title=None, title_loc=None, suptitle=None, ylabel=None, xlabel=None, add_empty_line_before_title=None, legend_title=None, legend=None):
    ax=axes.imshow(background_image, cmap='gray', aspect="auto")
    ax.axes.get_xaxis().set_visible(False)
    axes.grid(False)
    axes.set_yticklabels([])
    if ylabel:
        axes.set_ylabel(ylabel, fontsize=25, labelpad=5) # labelpad=5
    else:
        axes.set_ylabel('', fontsize=0.0001, labelpad=5) # labelpad=5
        
    if title_prefix:
        plot_title_=r"$\bf{" + title_prefix + "}$"
        print(plot_title_)
        if suptitle:
            plot_title_+=suptitle
            if plot_title:
                plot_title=plot_title_+'\n'+plot_title
            else:
                plot_title=plot_title_
        else:
            if plot_title:
                plot_title=plot_title_+plot_title
            else:
                plot_title=plot_title_
    if add_empty_line_before_title:
        if plot_title:
            plot_title='\n'+plot_title
        else:
            plot_title='\n'
    axes.set_title(plot_title, loc=title_loc, fontsize=25) # pad=3
    ax = sns.scatterplot(ax=axes, data=scatter_plot_data, x=x, y=y, hue=celltype_column_name, palette=palette, s=50)
    ax.invert_yaxis()
    if legend:
        sns.move_legend(
            ax,
            "lower center",
            bbox_to_anchor=[1.090, -0.400],
            ncol=3,
            title=legend_title,
            frameon=True,
            markerscale=3,
            fontsize=20,
        )
    else:
        ax.legend([],[], frameon=False)
    return ylabel
    
    