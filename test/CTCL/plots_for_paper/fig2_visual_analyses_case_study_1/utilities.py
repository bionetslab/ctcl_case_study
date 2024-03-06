import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import seaborn_image as isns

def _filter_dataframe_by_column_value_(dataframe, column_name, column_value):
    dataframe=dataframe[dataframe[column_name]==column_value]
    return dataframe

def _filter_dataframe_iteratively_by_list_of_column_values_(dataframe, column_name, list_of_column_values):
    filtered_results=[]
    for i in list_of_column_values: # or, <list_of_patient_ids>
        filtered_results.append(dataframe[dataframe[column_name]==i])
    filtered_results=pd.concat(filtered_results, axis=0)
    return filtered_results

def _generate_scatter_plot_with_overlay_(background_image, scatter_plot_data, x, y, celltype_column_name, title, savefig_name, palette='Set3', ncol=3):
    fig, axes = plt.subplots(figsize=(10,10))
    ax=axes.imshow(background_image, cmap='gray')
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax = sns.scatterplot(data=scatter_plot_data, x=x, y=y, hue=celltype_column_name, palette=palette, s=20)
    ax.set(xlabel=None, ylabel=None, title=title)
    ax.invert_yaxis()
    sns.move_legend(
        ax,
        "lower center",
        # bbox_to_anchor=(0.00, 0.00),
        ncol=3,
        title=None,
        frameon=True,
    )
    ax.grid(False) # Turns off grid.
    plt.savefig(savefig_name, format='svg', bbox_inches='tight')
    plt.show()
    plt.close()

def _generate_scatter_subplot_in_mosaic_(background_image, scatter_plot_data, x, y, celltype_column_name, subplot_axis_id, title, savefig_name, palette='Set3', ncol=3):
    # ax=isns.imgplot(background_image, cmap='gray', ax=subplot_axis_id)
    ax=subplot_axis_id.imshow(background_image, cmap='gray')
    # ax. margins(x=0)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax = sns.scatterplot(ax=subplot_axis_id, data=scatter_plot_data, x=x, y=y, hue=celltype_column_name, palette=palette, s=20)
    ax.set(xlabel=None, ylabel=None, title=title)
    ax.invert_yaxis()
    sns.move_legend(
        ax,
        "lower center",
        # bbox_to_anchor=(0.00, 0.00),
        ncol=3,
        title=None,
        frameon=True,
    )
    ax.grid(False) # Turns off grid.

def _generate_scatter_subplot_(background_image, scatter_plot_data, x, y, celltype_column_name, subplot_axis_id, title, savefig_name, palette='Set3', ncol=3):
    # ax=isns.imgplot(background_image, cmap='gray', ax=subplot_axis_id)
    ax=subplot_axis_id.imshow(background_image, cmap='gray')
    # ax. margins(x=0)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax = sns.scatterplot(ax=subplot_axis_id, data=scatter_plot_data, x=x, y=y, hue=celltype_column_name, palette=palette, s=10)
    ax.set(xlabel=None, ylabel=None)
    ax.set_title(title, fontsize=15)
    ax.invert_yaxis()
    sns.move_legend(
        ax,
        "lower center",
        # bbox_to_anchor=(0.00, 0.00),
        ncol=3,
        title=None,
        frameon=True,
    )
    ax.grid(False) # Turns off grid.
























