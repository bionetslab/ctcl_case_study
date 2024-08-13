import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
sample_wise_cell_types=pd.read_csv('result/sample-wise celltypes (HPA-based clustering).csv')
sample_wise_cell_types=sample_wise_cell_types.rename(columns={'celltype.1':'celltype'})
celltype_labels=sample_wise_cell_types.celltype
sample_wise_cell_types=sample_wise_cell_types.drop(columns=['x', 'y', 'sample_id',
'condition', 'celltype', 'cell_idx', 'cell_id', '_index_',
'condition.1', 'x.1', 'y.1'])
pca = PCA(n_components=3)
pca_result = pca.fit_transform(sample_wise_cell_types.values)
data_pca=pd.DataFrame(pca_result)
data_pca['celltype']=list(celltype_labels.iloc[:,0])
data_pca.to_csv('pca_projections.csv')
# import seaborn as sns
# sns.set_theme(style="whitegrid")
# import matplotlib.pyplot as plt
# ax = sns.scatterplot(
#     data=data_tsne, x=0, y=1, hue='species'
# )
# # plt.subplots_adjust(wspace=1.00, hspace=0.3)
# plt.savefig('t-SNE.pdf', format='pdf', bbox_inches='tight')
# plt.show()



