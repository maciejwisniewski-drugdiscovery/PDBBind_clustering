import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt

#  Outputs
clustalo_distance_matrix_heatmap_filepath = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/clustalo_matrices/clustalo_PDBBind_protein_sequences_distance_heatmap.png'

#  Load DataFrame with TYPES
dataframe_filepath='/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/dataframes/LP_PDBBind.csv'
dataframe = pd.read_csv(dataframe_filepath)
dataframe = dataframe.sort_values(by=['type']).reset_index()

#  Load Distance Matrix Tensor

if torch.cuda.is_available():
    print('GPU available')
    torch_distance_tensor_fp = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/clustalo_matrices/gpu_clustalo_PDBBind_protein_sequences_distance_matrix.pt'
else:
    print('GPU not available')
    torch_distance_tensor_fp = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/clustalo_matrices/clustalo_PDBBind_protein_sequences_distance_matrix.pt'

distance_tensor = torch.load(torch_distance_tensor_fp)
distance_tensor = distance_tensor.cpu()

#  Types of proteins
x_categories = dataframe['type'].unique()
y_categories = dataframe['type'].unique()
x_divisions = [0] + list(dataframe.groupby('type').size().cumsum())
y_divisions = x_divisions

#  HeatMap Creation
print('Heat Map Creation')
plt.figure()
plt.imshow(distance_tensor, cmap='hot', interpolation='nearest')

plt.xticks(np.array(x_divisions[:-1]) + np.diff(x_divisions), x_categories, rotation=45)
plt.yticks(np.array(y_divisions[:-1]) + np.diff(y_divisions), y_categories)

plt.title('Distance Matrix of Protein Sequences in PDBBind (Clustal Omega)')
plt.xlabel('Proteins')
plt.ylabel('Proteins')
plt.colorbar()

plt.savefig(clustalo_distance_matrix_heatmap_filepath)


