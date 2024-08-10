import seaborn as sns
from matplotlib.lines import Line2D

def _show_gray_image_(image, axes, title_prefix=None, plot_title=None, title_loc=None, xlabel=None, cmap=None, alpha=None, legend=None, legend_labels=None):
    ax=axes.imshow(image, cmap=cmap, aspect="auto", alpha=alpha)
    ax.axes.get_xaxis().set_visible(False)
    axes.grid(False)
    axes.set_yticklabels([])
    sns.set(font_scale = 1.2)
    sns.set_style("white")
    if title_prefix:
        if plot_title:
            plot_title=title_prefix + plot_title
        else:
            plot_title=title_prefix
    axes.set_title(plot_title, loc=title_loc, fontsize=35, fontweight='bold') # fontsize=25, pad=20
    axes.grid(False)
    
    if legend:
        # ========== Generate legend: ===========
        legend_elements = [Line2D([0], [0], color='gray', lw=10, label=legend_labels[0]),
                           Line2D([0], [0], color='red', lw=10, label=legend_labels[1])]
        # axes.legend(handles=legend_elements, bbox_to_anchor=(1.00, 0.90), loc='lower left', fontsize=23)
        axes.legend(handles=legend_elements, loc='upper right', fontsize=30)
    axes.grid(False)
    