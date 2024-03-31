import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import seaborn_image as isns

def _generate_scatter_subplot_(background_image, scatter_plot_data, x, y, celltype_column_name, axes, palette='Set3', ncol=3, title=None):
    # # ax=isns.imgplot(background_image, cmap='gray', ax=subplot_axis_id)
    # ax=subplot_axis_id.imshow(background_image, cmap='gray')
    # # ax. margins(x=0)
    # ax.axes.get_xaxis().set_visible(False)
    # ax.axes.get_yaxis().set_visible(False)
    if title_prefix:
        plot_title_=r"$\bf{" + title_prefix + "}$"
        if suptitle:
            plot_title_+=suptitle
            plot_title=plot_title_+'\n'+plot_title
        else:
            plot_title=plot_title_+plot_title
    ax.set_title(plot_title, loc=title_loc, pad=2)
    ax = sns.scatterplot(ax=axes, data=scatter_plot_data, x=x, y=y, hue=celltype_column_name, palette=palette, s=100)
    ax.set(xlabel=None, ylabel=None)
    ax.set_title(title, fontsize=15)
    ax.invert_yaxis()
    # subplot_axis_id.imshow(background_image, cmap='gray')
    sns.move_legend(
        ax,
        "lower center",
        # bbox_to_anchor=(0.00, 0.00),
        ncol=3,
        title=None,
        frameon=True,
        markerscale=3,
        fontsize=20,
    )
    ax.grid(False) # Turns off grid.
    # ax.set_aspect(0.25)
    ax.axis('off')
    # ax.margins(y=0)
    # ax=subplot_axis_id.imshow(background_image, cmap='gray')
    # # ax.legend.set_title(None)
    # legend = axes.legend(frameon=True)
    # for legend_handle in legend.legendHandles:
    #     legend_handle._legmarker.set_markersize(9)
    
    # plt.legend(title=None, fontsize='50', title_fontsize='14')