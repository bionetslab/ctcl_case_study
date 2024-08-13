from sklearn.manifold import TSNE
import pandas as pd
sample_wise_cell_types=pd.read_csv('result/sample-wise celltypes (HPA-based clustering).csv')
sample_wise_cell_types=sample_wise_cell_types.rename(columns={'celltype.1':'celltype'})
celltype_labels=sample_wise_cell_types.celltype
sample_wise_cell_types=sample_wise_cell_types.drop(columns=['x', 'y', 'sample_id',
'condition', 'celltype', 'cell_idx', 'cell_id', '_index_',
'condition.1', 'x.1', 'y.1'])
tsne = TSNE(n_components=2, random_state=0)
projections = tsne.fit_transform(sample_wise_cell_types)
data_tsne=pd.DataFrame(projections)
data_tsne['celltype']=list(celltype_labels.iloc[:,0])
data_tsne.to_csv('tsne_projections.csv')
# import seaborn as sns
# sns.set_theme(style="whitegrid")
# import matplotlib.pyplot as plt
# ax = sns.scatterplot(
#     data=data_tsne, x=0, y=1, hue='species'
# )
# # plt.subplots_adjust(wspace=1.00, hspace=0.3)
# plt.savefig('t-SNE.pdf', format='pdf', bbox_inches='tight')
# plt.show()



