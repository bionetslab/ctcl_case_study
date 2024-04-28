import seaborn as sns
import matplotlib.pyplot as plt

def plot_heatmap(data, genes, figsize=(9,6), outfile=None, title=None):
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(data=data[genes],ax=ax)
    ax.set(ylabel='')
    if not title is None:
        ax.set_title(title)
    fig.tight_layout()
    if not outfile is None:
        fig.savefig(outfile)